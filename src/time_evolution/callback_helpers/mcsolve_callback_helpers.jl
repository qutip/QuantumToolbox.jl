#=
Helper functions for the mcsolve callbacks.
=#

struct SaveFuncMCSolve{TE,IT,TEXPV}
    e_ops::TE
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncMCSolve)(integrator) = _save_func_mcsolve(integrator, f.e_ops, f.iter, f.expvals)

struct LindbladJump{T1,T2}
    c_ops::T1
    c_ops_herm::T2
end

(f::LindbladJump)(integrator) = _lindblad_jump_affect!(integrator, f.c_ops, f.c_ops_herm)

##

function _save_func_mcsolve(integrator, e_ops, iter, expvals)
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
        # We are implicitly saying that we don't have a `ProgressBar`
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb1, kwargs.callback),)) :
            merge(kwargs, (callback = cb1,))
        return kwargs2
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(tlist))

        _save_affect! = SaveFuncMCSolve(get_data.(e_ops), Ref(1), expvals)
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

function _mcsolve_get_expvals(sol::AbstractODESolution)
    cb = NamedTuple(sol.prob.kwargs).callback
    if _mcsolve_has_discrete_callbacks(cb)
        return _mcsolve_get_expvals(cb)
    else
        return nothing
    end
end
function _mcsolve_get_expvals(cb::CallbackSet)
    idx = _mcsolve_has_continuous_jump(cb) ? 1 : 2
    _cb = cb.discrete_callbacks[idx]
    return _mcsolve_get_expvals(_cb)
end
_mcsolve_get_expvals(cb::DiscreteCallback) =
    if cb.affect! isa SaveFuncMCSolve
        return cb.affect!.expvals
    else
        nothing
    end
_mcsolve_get_expvals(cb::ContinuousCallback) = nothing

#=
    _mcsolve_callbacks_new_iter_expvals(prob, tlist)

Return the same callbacks of the `prob`, but with the `iter` variable reinitialized to 1 and the `expvals` variable reinitialized to a new matrix.
=#
function _mcsolve_callbacks_new_iter_expvals(prob, tlist)
    cb = prob.kwargs[:callback]
    return _mcsolve_callbacks_new_iter_expvals(cb, tlist)
end
function _mcsolve_callbacks_new_iter_expvals(cb::CallbackSet, tlist)
    cb_continuous = cb.continuous_callbacks
    cb_discrete = cb.discrete_callbacks

    if _mcsolve_has_continuous_jump(cb)
        idx = 1
        e_ops = cb_discrete[idx].affect!.e_ops
        expvals = similar(cb_discrete[idx].affect!.expvals)
        _save_affect! = SaveFuncMCSolve(e_ops, Ref(1), expvals)
        cb_save = PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
        return CallbackSet(cb_continuous..., cb_save, cb_discrete[2:end]...)
    else
        idx = 2
        e_ops = cb_discrete[idx].affect!.e_ops
        expvals = similar(cb_discrete[idx].affect!.expvals)
        _save_affect! = SaveFuncMCSolve(e_ops, Ref(1), expvals)
        cb_save = PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
        return CallbackSet(cb_continuous..., cb_discrete[1], cb_save, cb_discrete[3:end]...)
    end
end
_mcsolve_callbacks_new_iter_expvals(cb::ContinuousCallback, tlist) = cb # It is only the continuous LindbladJump callback  
_mcsolve_callbacks_new_iter_expvals(cb::DiscreteCallback, tlist) = cb # It is only the discrete LindbladJump callback

_mcsolve_has_discrete_callbacks(cb::CallbackSet) = length(cb.discrete_callbacks) > 0
_mcsolve_has_discrete_callbacks(cb::ContinuousCallback) = false
_mcsolve_has_discrete_callbacks(cb::DiscreteCallback) = true

_mcsolve_has_continuous_jump(cb::CallbackSet) =
    (length(cb.continuous_callbacks) > 0) && (cb.continuous_callbacks[1].affect! isa LindbladJump)
_mcsolve_has_continuous_jump(cb::ContinuousCallback) = true
_mcsolve_has_continuous_jump(cb::DiscreteCallback) = false
