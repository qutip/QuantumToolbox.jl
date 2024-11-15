#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

########## SESOLVE ##########

struct SaveFuncSESolve{T1,T2}
    e_ops::T1
    is_empty_e_ops::T2
end

(f::SaveFuncSESolve)(integrator) = _save_func_sesolve(integrator, f.e_ops, f.is_empty_e_ops)
(f::SaveFuncSESolve{Nothing})(integrator) = _save_func_sesolve(integrator)

# When e_ops is Nothing
function _save_func_sesolve(integrator)
    next!(integrator.p.progr)
    u_modified!(integrator, false)
    return nothing
end

# When e_ops is a list of operators
function _save_func_sesolve(integrator, e_ops, is_empty_e_ops)
    expvals = integrator.p.expvals
    progr = integrator.p.progr
    if !is_empty_e_ops
        ψ = integrator.u
        _expect = op -> dot(ψ, get_data(op), ψ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    return _save_func_sesolve(integrator)
end

function _generate_sesolve_callback(e_ops, tlist)
    is_empty_e_ops = e_ops isa Nothing ? true : isempty(e_ops)
    _save_affect! = SaveFuncSESolve(e_ops, is_empty_e_ops)
    return PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
end

########## MCSOLVE ##########

struct SaveFuncMCSolve{T1,T2}
    e_ops::T1
    is_empty_e_ops::T2
end

(f::SaveFuncMCSolve)(integrator) = _save_func_mcsolve(integrator, f.e_ops, f.is_empty_e_ops)

struct LindbladJumpAffect!{T1,T2}
    c_ops::T1
    c_ops_herm::T2
end

(f::LindbladJumpAffect!)(integrator) = _lindblad_jump_affect!(integrator, f.c_ops, f.c_ops_herm)

##

function _save_func_mcsolve(integrator, e_ops, is_empty_e_ops)
    expvals = integrator.p.expvals
    progr = integrator.p.progr
    cache_mc = integrator.p.mcsolve_params.cache_mc
    if !is_empty_e_ops
        copyto!(cache_mc, integrator.u)
        normalize!(cache_mc)
        ψ = cache_mc
        _expect = op -> dot(ψ, get_data(op), ψ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    next!(progr)
    u_modified!(integrator, false)
    return nothing
end

function _generate_mcsolve_kwargs(e_ops, tlist, c_ops, jump_callback, kwargs)
    c_ops_data = get_data.(c_ops)
    c_ops_herm_data = map(op -> op' * op, c_ops_data)

    _affect! = LindbladJumpAffect!(c_ops_data, c_ops_herm_data)

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
        is_empty_e_ops = isempty(e_ops)
        _save_affect! = SaveFuncMCSolve(e_ops, is_empty_e_ops)
        cb2 = PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
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
