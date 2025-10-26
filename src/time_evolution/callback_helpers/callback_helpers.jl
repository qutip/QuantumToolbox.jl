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

    # Ensure that the noise is stored in tlist. # TODO: Fix this directly in DiffEqNoiseProcess.jl
    # See https://github.com/SciML/DiffEqNoiseProcess.jl/issues/214 for example
    tstops = haskey(kwargs, :tstops) ? unique!(sort!(vcat(tlist, kwargs.tstops))) : tlist
    kwargs2 = merge(kwargs, (tstops = tstops,))

    if SF === SaveFuncSSESolve
        cb_normalize = _ssesolve_generate_normalize_cb()
        return _merge_kwargs_with_callback(kwargs2, CallbackSet(cb_normalize, cb_save))
    end

    return _merge_kwargs_with_callback(kwargs2, cb_save)
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

    progr =
        getVal(progress_bar) ?
        Progress(length(tlist), showspeed = true, enabled = getVal(progress_bar), desc = _get_progress_desc(method)) :
        nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))

    _save_func = method(e_ops_data, progr, Ref(1), expvals)
    return FunctionCallingCallback(_save_func, funcat = tlist)
end

function _generate_stochastic_save_callback(e_ops, sc_ops, tlist, store_measurement, progress_bar, method)
    e_ops_data = e_ops isa Nothing ? nothing : _get_e_ops_data(e_ops, method)
    m_ops_data = _get_m_ops_data(sc_ops, method)

    progr = getVal(progress_bar) ? Progress(length(tlist), showspeed = true, enabled = getVal(progress_bar)) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))
    m_expvals = getVal(store_measurement) ? Array{Float64}(undef, length(sc_ops), length(tlist) - 1) : nothing

    _save_func_cache = Array{Float64}(undef, length(sc_ops))
    _save_func =
        method(store_measurement, e_ops_data, m_ops_data, progr, Ref(1), expvals, m_expvals, tlist, _save_func_cache)
    return FunctionCallingCallback(_save_func, funcat = tlist)
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
        return cb.affect!.func.m_expvals
    end
end

#=
    With this function we extract the e_ops from the SaveFuncMCSolve `affect!` function of the callback of the integrator.
    This callback can only be a FunctionCallingCallback (DiscreteCallback).
=#
function _get_e_ops(integrator::AbstractODEIntegrator, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cb = _get_save_callback(integrator, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.func.e_ops
    end
end

# Get the e_ops from a given AbstractODESolution. Valid for `sesolve`, `mesolve` and `ssesolve`.
function _get_expvals(sol::AbstractODESolution, method::Type{SF}) where {SF<:AbstractSaveFunc}
    cb = _get_save_callback(sol, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.func.expvals
    end
end

#=
    _get_save_callback

Return the Callback that is responsible for saving the expectation values of the system.
=#
function _get_save_callback(sol::AbstractODESolution, method::Type{SF}) where {SF<:AbstractSaveFunc}
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    if hasproperty(kwargs, :callback) && !isnothing(kwargs.callback)
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
    if typeof(cb.affect!) <: FunctionCallingAffect && typeof(cb.affect!.func) <: AbstractSaveFunc
        return cb
    end
    return nothing
end
_get_save_callback(cb::ContinuousCallback, ::Type{SF}) where {SF<:AbstractSaveFunc} = nothing

_get_save_callback_idx(cb, method) = 1

# %% ------------ Noise Measurement Helpers ------------ %%

# TODO: To improve. See https://github.com/SciML/DiffEqNoiseProcess.jl/issues/214
function _homodyne_dWdt!(dWdt_cache, integrator, tlist, iter)
    idx = findfirst(>=(tlist[iter[]-1]), integrator.W.t)

    # We are assuming that the last element is tlist[iter[]]
    @inbounds dWdt_cache .= (integrator.W.u[end] .- integrator.W.u[idx]) ./ (integrator.W.t[end] - integrator.W.t[idx])

    return nothing
end
