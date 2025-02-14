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

function _generate_stochastic_kwargs(
    e_ops,
    sc_ops,
    progress_bar,
    tlist,
    store_measurement,
    kwargs,
    method::Type{SF},
) where {SF<:AbstractSaveFunc}
    cb_save = _generate_stochastic_save_callback(e_ops, sc_ops, tlist, store_measurement, progress_bar, method)

    if SF === SaveFuncSSESolve
        cb_normalize = _ssesolve_generate_normalize_cb()
        return _merge_kwargs_with_callback(kwargs, CallbackSet(cb_normalize, cb_save))
    end

    return _merge_kwargs_with_callback(kwargs, cb_save)
end
_generate_stochastic_kwargs(
    e_ops::Nothing,
    sc_ops,
    progress_bar::Val{false},
    tlist,
    store_measurement::Val{false},
    kwargs,
    method::Type{SF},
) where {SF<:AbstractSaveFunc} = kwargs

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

function _generate_stochastic_save_callback(e_ops, sc_ops, tlist, store_measurement, progress_bar, method)
    e_ops_data = e_ops isa Nothing ? nothing : _get_e_ops_data(e_ops, method)
    m_ops_data = _get_m_ops_data(sc_ops, method)

    progr = getVal(progress_bar) ? ProgressBar(length(tlist), enable = getVal(progress_bar)) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))
    m_expvals = getVal(store_measurement) ? Array{Float64}(undef, length(sc_ops), length(tlist) - 1) : nothing

    _save_affect! = method(store_measurement, e_ops_data, m_ops_data, progr, Ref(1), expvals, m_expvals)
    return PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
end

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
    To extract the measurement outcomes of a stochastic solver
=#
function _get_m_expvals(integrator::AbstractODESolution, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cb = _get_save_callback(integrator, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.m_expvals
    end
end

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

# %% ------------ Noise Measurement Helpers ------------ %%

# TODO: Add some cache mechanism to avoid memory allocations
function _homodyne_dWdt(integrator)
    @inbounds _dWdt = (integrator.W.u[end] .- integrator.W.u[end-1]) ./ (integrator.W.t[end] - integrator.W.t[end-1])

    return _dWdt
end
