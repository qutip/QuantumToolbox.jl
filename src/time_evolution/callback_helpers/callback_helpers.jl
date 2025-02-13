#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

##

abstract type AbstractSaveFunc end

# Multiple dispatch depending on the progress_bar and e_ops types
function _generate_se_me_kwargs(e_ops, progress_bar, tlist, kwargs, method)
    cb = _generate_save_callback(e_ops, tlist, progress_bar, method)
    return _merge_kwargs_with_callback(kwargs, cb)
end
_generate_se_me_kwargs(e_ops::Nothing, progress_bar::Val{false}, tlist, kwargs, method) = kwargs

function _merge_kwargs_with_callback(kwargs, cb)
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb, kwargs.callback),)) :
        merge(kwargs, (callback = cb,))

    return kwargs2
end

function _generate_save_callback(e_ops, tlist, progress_bar, method)
    e_ops_data = e_ops isa Nothing ? nothing : _get_e_ops_data(e_ops, method)

    progr = getVal(progress_bar) ? ProgressBar(length(tlist), enable = getVal(progress_bar)) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))

    _save_affect! = method(e_ops_data, progr, Ref(1), expvals)
    return PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
end

_generate_mesolve_e_op(op) = mat2vec(adjoint(get_data(op)))

##

# When e_ops is Nothing. Common for all solvers
function _save_func(integrator, progr)
    next!(progr)
    u_modified!(integrator, false)
    return nothing
end

# When progr is Nothing. Common for all solvers
function _save_func(integrator, progr::Nothing)
    u_modified!(integrator, false)
    return nothing
end

##

#=
    With this function we extract the e_ops from the SaveFuncMCSolve `affect!` function of the callback of the integrator.
    This callback can only be a PresetTimeCallback (DiscreteCallback).
=#
function _get_e_ops(integrator::AbstractODEIntegrator, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cb = _get_save_callback(integrator, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.e_ops
    end
end

# Get the e_ops from a given AbstractODESolution. Valid for `sesolve`, `mesolve` and `ssesolve`.
function _get_expvals(sol::AbstractODESolution, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cb = _get_save_callback(sol, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.expvals
    end
end

#=
    _get_save_callback

Return the Callback that is responsible for saving the expectation values of the system.
=#
function _get_save_callback(sol::AbstractODESolution, method::Type{SF}) where {SF<:AbstractSaveFunc}
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    if hasproperty(kwargs, :callback)
        return _get_save_callback(kwargs.callback, method)
    else
        return nothing
    end
end
_get_save_callback(integrator::AbstractODEIntegrator, method::Type{SF}) where {SF<:AbstractSaveFunc} =
    _get_save_callback(integrator.opts.callback, method)
function _get_save_callback(cb::CallbackSet, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cbs_discrete = cb.discrete_callbacks
    if length(cbs_discrete) > 0
        idx = _get_save_callback_idx(cb, method)
        _cb = cb.discrete_callbacks[idx]
        return _get_save_callback(_cb, method)
    else
        return nothing
    end
end
function _get_save_callback(cb::DiscreteCallback, ::Type{SF}) where {SF<:AbstractSaveFunc}
    if typeof(cb.affect!) <: AbstractSaveFunc
        return cb
    end
    return nothing
end
_get_save_callback(cb::ContinuousCallback, ::Type{SF}) where {SF<:AbstractSaveFunc} = nothing

_get_save_callback_idx(cb, method) = 1
