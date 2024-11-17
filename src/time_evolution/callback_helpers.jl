#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

########## SESOLVE ##########

struct SaveFuncSESolve{TE,PT<:Union{Nothing,ProgressBar},IT}
    e_ops::TE
    progr::PT
    iter::IT
end

(f::SaveFuncSESolve)(integrator) = _save_func_sesolve(integrator, f.e_ops, f.progr, f.iter)
(f::SaveFuncSESolve{Nothing})(integrator) = _save_func_sesolve(integrator, f.progr)

##

# When e_ops is Nothing
function _save_func_sesolve(integrator, progr)
    next!(progr)
    u_modified!(integrator, false)
    return nothing
end

# When progr is Nothing
function _save_func_sesolve(integrator, progr::Nothing)
    u_modified!(integrator, false)
    return nothing
end

# When e_ops is a list of operators
function _save_func_sesolve(integrator, e_ops, progr, iter)
    expvals = integrator.p.expvals
    ψ = integrator.u
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func_sesolve(integrator, progr)
end

function _generate_sesolve_callback(e_ops, tlist, progress_bar)
    e_ops_data = e_ops isa Nothing ? nothing : get_data.(e_ops)

    progr = getVal(progress_bar) ? ProgressBar(length(tlist), enable = getVal(progress_bar)) : nothing

    _save_affect! = SaveFuncSESolve(e_ops_data, progr, Ref(1))
    return PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
end

########## MCSOLVE ##########

struct SaveFuncMCSolve{TE,IT}
    e_ops::TE
    iter::IT
end

(f::SaveFuncMCSolve)(integrator) = _save_func_mcsolve(integrator, f.e_ops, f.iter)

struct LindbladJump{T1,T2}
    c_ops::T1
    c_ops_herm::T2
end

(f::LindbladJump)(integrator) = _lindblad_jump_affect!(integrator, f.c_ops, f.c_ops_herm)

##

function _save_func_mcsolve(integrator, e_ops, iter)
    expvals = integrator.p.expvals
    cache_mc = integrator.p.mcsolve_params.cache_mc

    copyto!(cache_mc, integrator.u)
    normalize!(cache_mc)
    ψ = cache_mc
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    u_modified!(integrator, false)
    return nothing
end

function _generate_mcsolve_kwargs(e_ops, tlist, c_ops, jump_callback, kwargs)
    c_ops_data = get_data.(c_ops)
    c_ops_herm_data = map(op -> op' * op, c_ops_data)

    _affect! = LindbladJump(c_ops_data, c_ops_herm_data)

    if jump_callback isa DiscreteLindbladJumpCallback
        cb1 = DiscreteCallback(_mcsolve_discrete_condition, _affect!, save_positions = (false, false))
    else
        cb1 = ContinuousCallback(
            _mcsolve_continuous_condition,
            _affect!,
            nothing,
            interp_points = jump_callback.interp_points,
            save_positions = (false, false),
        )
    end

    if e_ops isa Nothing
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb1, kwargs.callback),)) :
            merge(kwargs, (callback = cb1,))
        return kwargs2
    else
        _save_affect! = SaveFuncMCSolve(get_data.(e_ops), Ref(1))
        cb2 = _PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb1, cb2, kwargs.callback),)) :
            merge(kwargs, (callback = CallbackSet(cb1, cb2),))
        return kwargs2
    end
end

function _lindblad_jump_affect!(integrator, c_ops, c_ops_herm)
    params = integrator.p
    cache_mc = params.mcsolve_params.cache_mc
    weights_mc = params.mcsolve_params.weights_mc
    cumsum_weights_mc = params.mcsolve_params.cumsum_weights_mc
    random_n = params.mcsolve_params.random_n
    jump_times = params.mcsolve_params.jump_times
    jump_which = params.mcsolve_params.jump_which
    jump_times_which_idx = params.mcsolve_params.jump_times_which_idx
    traj_rng = params.mcsolve_params.traj_rng
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        weights_mc[i] = real(dot(ψ, c_ops_herm[i], ψ))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    r = rand(traj_rng) * sum(real, weights_mc)
    collapse_idx = getindex(1:length(weights_mc), findfirst(x -> real(x) > r, cumsum_weights_mc))
    mul!(cache_mc, c_ops[collapse_idx], ψ)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)

    @inbounds random_n[1] = rand(traj_rng)

    @inbounds idx = round(Int, real(jump_times_which_idx[1]))
    @inbounds jump_times[idx] = integrator.t
    @inbounds jump_which[idx] = collapse_idx
    @inbounds jump_times_which_idx[1] += 1
    @inbounds if real(jump_times_which_idx[1]) > length(jump_times)
        resize!(jump_times, length(jump_times) + JUMP_TIMES_WHICH_INIT_SIZE)
        resize!(jump_which, length(jump_which) + JUMP_TIMES_WHICH_INIT_SIZE)
    end
    u_modified!(integrator, true)
    return nothing
end

_mcsolve_continuous_condition(u, t, integrator) =
    @inbounds real(integrator.p.mcsolve_params.random_n[1]) - real(dot(u, u))

_mcsolve_discrete_condition(u, t, integrator) =
    @inbounds real(dot(u, u)) < real(integrator.p.mcsolve_params.random_n[1])

#=
With this function we extract the c_ops and c_ops_herm from the LindbladJump `affect!` function of the callback of the integrator.
This callback can be a DiscreteLindbladJumpCallback or a ContinuousLindbladJumpCallback.
=#
function _mcsolve_get_c_ops(integrator::AbstractODEIntegrator)
    cb_set = integrator.opts.callback # This is supposed to be a CallbackSet
    (cb_set isa CallbackSet) || throw(ArgumentError("The callback must be a CallbackSet."))
    cb = isempty(cb_set.continuous_callbacks) ? cb_set.discrete_callback[1] : cb_set.continuous_callbacks[1]
    return cb.affect!.c_ops, cb.affect!.c_ops_herm
end

#=
With this function we extract the e_ops from the SaveFuncMCSolve `affect!` function of the callback of the integrator.
This callback can only be a PresetTimeCallback (DiscreteCallback).
=#
function _mcsolve_get_e_ops(integrator::AbstractODEIntegrator)
    cb_set = integrator.opts.callback # This is supposed to be a CallbackSet
    (cb_set isa CallbackSet) || throw(ArgumentError("The callback must be a CallbackSet."))
    cb = length(cb_set.continuous_callbacks) > 0 ? cb_set.discrete_callbacks[1] : cb_set.discrete_callbacks[2]
    return cb.affect!.e_ops
end

#=
    _mcsolve_callbacks_new_iter(prob, tlist)

Return the same callbacks of the `prob`, but with the `iter` variable reinitialized to 1.
=#
function _mcsolve_callbacks_new_iter(prob, tlist)
    cb = prob.kwargs[:callback]
    return _mcsolve_callbacks_new_iter(cb, tlist)
end
function _mcsolve_callbacks_new_iter(cb::CallbackSet, tlist)
    cb_continuous = cb.continuous_callbacks
    cb_discrete = cb.discrete_callbacks

    if length(cb_continuous) > 0
        idx = 1
        e_ops = cb_discrete[idx].affect!.e_ops
        _save_affect! = SaveFuncMCSolve(e_ops, Ref(1))
        cb_save = _PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
        return CallbackSet(cb_continuous..., cb_save, cb_discrete[2:end]...)
    else
        idx = 2
        e_ops = cb_discrete[idx].affect!.e_ops
        _save_affect! = SaveFuncMCSolve(e_ops, Ref(1))
        cb_save = _PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
        return CallbackSet(cb_continuous..., cb_discrete[1], cb_save, cb_discrete[3:end]...)
    end
end
_mcsolve_callbacks_new_iter(cb::ContinuousCallback, tlist) = cb
_mcsolve_callbacks_new_iter(cb::DiscreteCallback, tlist) = cb

## Temporary function to avoid errors. Waiting for the PR In DiffEqCallbacks.jl to be merged.

import SciMLBase: INITIALIZE_DEFAULT, add_tstop!

function _PresetTimeCallback(
    tstops,
    user_affect!;
    initialize = INITIALIZE_DEFAULT,
    filter_tstops = true,
    sort_inplace = false,
    kwargs...,
)
    if !(tstops isa AbstractVector) && !(tstops isa Number)
        throw(ArgumentError("tstops must either be a number or a Vector. Was $tstops"))
    end

    tstops = tstops isa Number ? [tstops] : (sort_inplace ? sort!(tstops) : sort(tstops))

    condition = let
        function (u, t, integrator)
            if hasproperty(integrator, :dt)
                insorted(t, tstops) && (integrator.t - integrator.dt) != integrator.t
            else
                insorted(t, tstops)
            end
        end
    end

    # Initialization: first call to `f` should be *before* any time steps have been taken:
    initialize_preset = function (c, u, t, integrator)
        initialize(c, u, t, integrator)

        if filter_tstops
            tdir = integrator.tdir
            tspan = integrator.sol.prob.tspan
            _tstops = tstops[@. tdir * tspan[1] < tdir * tstops < tdir * tspan[2]]
        else
            _tstops = tstops
        end
        for tstop in _tstops
            add_tstop!(integrator, tstop)
        end
        if t ∈ tstops
            user_affect!(integrator)
        end
    end
    return DiscreteCallback(condition, user_affect!; initialize = initialize_preset, kwargs...)
end
