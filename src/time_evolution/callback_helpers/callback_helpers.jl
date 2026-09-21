#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

##

abstract type AbstractSaveFunc end

function _merge_tstops(kwargs, skip_tstops::Bool, tlist)
    if skip_tstops
        return kwargs
    else
        tstops = haskey(kwargs, :tstops) ? unique!(sort!(vcat(tlist, kwargs.tstops))) : tlist
        return merge(kwargs, (tstops = tstops,))
    end
end

# Multiple dispatch depending on the progress_bar and e_ops types
function _generate_se_me_kwargs(e_ops, progress_bar, tlist, kwargs, method, ::Type{T}) where {T <: Number}
    cb = _generate_save_callback(e_ops, tlist, progress_bar, method, T)
    return _merge_kwargs_with_callback(kwargs, cb)
end
_generate_se_me_kwargs(e_ops::Nothing, progress_bar::Val{false}, tlist, kwargs, method, ::Type{T}) where {T <: Number} =
    kwargs

function _generate_stochastic_kwargs(
        e_ops,
        sc_ops,
        progress_bar,
        tlist,
        store_measurement,
        kwargs,
        method::Type{SF},
        ::Type{T},
    ) where {SF <: AbstractSaveFunc, T <: Number}
    cb_save = _generate_stochastic_save_callback(e_ops, sc_ops, tlist, store_measurement, progress_bar, method, T)

    # Ensure that the noise is stored in tlist. # TODO: Fix this directly in DiffEqNoiseProcess.jl
    # See https://github.com/SciML/DiffEqNoiseProcess.jl/issues/214 for example
    kwargs2 = _merge_tstops(kwargs, !getVal(store_measurement), tlist)

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
    ::Type{T},
) where {SF <: AbstractSaveFunc, T <: Number} = kwargs

function _merge_kwargs_with_callback(kwargs, cb)
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb, kwargs.callback),)) :
        merge(kwargs, (callback = cb,))

    return kwargs2
end

function _generate_save_callback(e_ops, tlist, progress_bar, method, ::Type{T}) where {T <: Number}
    e_ops_data = e_ops isa Nothing ? nothing : _get_e_ops_data(e_ops, method)

    progr =
        getVal(progress_bar) ?
        Progress(
            length(tlist);
            enabled = getVal(progress_bar),
            desc = _get_progress_desc(method),
            settings.ProgressMeterKWARGS...,
        ) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{_complex_float_type(T)}(undef, length(e_ops), length(tlist))

    _save_func = method(e_ops_data, progr, Ref(1), expvals)
    return FunctionCallingCallback(_save_func, funcat = tlist)
end

function _generate_stochastic_save_callback(
        e_ops,
        sc_ops,
        tlist,
        store_measurement,
        progress_bar,
        method,
        ::Type{T},
    ) where {T <: Number}
    e_ops_data = e_ops isa Nothing ? nothing : _get_e_ops_data(e_ops, method)
    m_ops_data = _get_m_ops_data(sc_ops, method)

    progr =
        getVal(progress_bar) ?
        Progress(length(tlist); enabled = getVal(progress_bar), settings.ProgressMeterKWARGS...) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{_complex_float_type(T)}(undef, length(e_ops), length(tlist))
    m_expvals = getVal(store_measurement) ? Array{_float_type(T)}(undef, length(sc_ops), length(tlist) - 1) : nothing

    _save_func_cache = Array{_float_type(T)}(undef, length(sc_ops))
    _save_func =
        method(store_measurement, e_ops_data, m_ops_data, progr, Ref(1), expvals, m_expvals, tlist, _save_func_cache)
    return FunctionCallingCallback(_save_func, funcat = tlist)
end

##

# When e_ops is Nothing. Common for all solvers
function _save_func(integrator, progr)
    next!(progr)
    derivative_discontinuity!(integrator, false)
    return nothing
end

# When progr is Nothing. Common for all solvers
function _save_func(integrator, progr::Nothing)
    derivative_discontinuity!(integrator, false)
    return nothing
end

##

#=
    With this function we extract the e_ops from the SaveFuncMCSolve `affect!` function of the callback of the integrator.
    This callback can only be a FunctionCallingCallback (DiscreteCallback).
=#
function _get_e_ops(integrator::AbstractODEIntegrator, method::Type{SF}) where {SF <: AbstractSaveFunc}
    cb = _get_save_callback(integrator, method)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.func.e_ops
    end
end

# Get the e_ops from a given AbstractODESolution. Valid for `sesolve`, `mesolve` and `mcsolve`.
#
# Note: depending on the solver internals (e.g. the SDE/ODE integrator), the callback stored in `sol`
# may have had its concrete type erased (e.g. `discrete_callbacks` stored as a `Vector{Any}` instead
# of a `Tuple`) to reduce compilation. This makes `cb.affect!.func.expvals` inferred as `Any`. Since we
# know the element type of `expvals` must match the (complex) element type of the solution itself, we
# recover type stability with an explicit type assertion.
function _get_expvals(sol::AbstractODESolution, method::Type{SF}) where {SF <: AbstractSaveFunc}
    cb = _get_save_callback(sol, method)
    return _get_expvals(cb, _complex_float_type(eltype(eltype(sol.u))))
end
_get_expvals(cb::Nothing, ::Type{CT}) where {CT <: Number} = nothing
_get_expvals(cb, ::Type{CT}) where {CT <: Number} = cb.affect!.func.expvals::Union{Nothing, Matrix{CT}}

#=
    Variants of `_get_expvals`/`_get_m_expvals` for the ensemble-based stochastic solvers
    (`ssesolve`/`smesolve`), where the caller already knows at compile time (from the type of
    `e_ops`/`store_measurement`, propagated as a `Val`) whether the result must be `nothing` or a
    concrete `Matrix`. Dispatching on that `Val`, instead of on a runtime `isnothing` check on the
    (possibly type-erased) callback, avoids `Union{Nothing, ...}` creeping into the return type of
    `ssesolve`/`smesolve` and breaking `@inferred`.
=#
_get_expvals(sol::AbstractODESolution, method::Type{SF}, ::Val{false}) where {SF <: AbstractSaveFunc} = nothing
function _get_expvals(sol::AbstractODESolution, method::Type{SF}, ::Val{true}) where {SF <: AbstractSaveFunc}
    cb = _get_save_callback(sol, method)
    return cb.affect!.func.expvals::Matrix{_complex_float_type(eltype(eltype(sol.u)))}
end

_get_m_expvals(sol::AbstractODESolution, method::Type{SF}, ::Val{false}) where {SF <: AbstractSaveFunc} = nothing
function _get_m_expvals(sol::AbstractODESolution, method::Type{SF}, ::Val{true}) where {SF <: AbstractSaveFunc}
    cb = _get_save_callback(sol, method)
    return cb.affect!.func.m_expvals::Matrix{_float_type(eltype(eltype(sol.u)))}
end

#=
    Stack per-trajectory `expvals`/`m_expvals` across an ensemble solution, given a compile-time `Val`
    flag (known from `e_ops`/`store_measurement` at the call site) stating whether they were requested.
    See the note above regarding why dispatching on `Val` (rather than a runtime check) is necessary.

    Note: these use an explicit loop into a preallocated `Vector`, rather than `map` with a closure
    over `method`. A closure capturing a `Type{SF} where {SF <: AbstractSaveFunc}`-typed argument is
    inferred as `Any` inside the closure (the argument's declared type is the abstract `where`-bound,
    not the concrete type it is instantiated with), which would make `_get_expvals`/`_get_m_expvals`
    dispatch dynamically and infer as `Any` again, undoing the fix above.
=#
function _stack_traj_expvals(::Val{true}, sol, method::Type{SF}) where {SF <: AbstractSaveFunc}
    CT = _complex_float_type(eltype(eltype(sol.u)))
    v = Vector{Matrix{CT}}(undef, length(sol.u))
    for i in eachindex(sol.u)
        v[i] = _get_expvals(sol.u[i], method, Val(true))
    end
    return stack(v, dims = 2)
end
_stack_traj_expvals(::Val{false}, sol, method::Type{SF}) where {SF <: AbstractSaveFunc} = nothing

function _stack_traj_m_expvals(::Val{true}, sol, method::Type{SF}) where {SF <: AbstractSaveFunc}
    FT = _float_type(eltype(eltype(sol.u)))
    v = Vector{Matrix{FT}}(undef, length(sol.u))
    for i in eachindex(sol.u)
        v[i] = _get_m_expvals(sol.u[i], method, Val(true))
    end
    return stack(v, dims = 2)
end
_stack_traj_m_expvals(::Val{false}, sol, method::Type{SF}) where {SF <: AbstractSaveFunc} = nothing

#=
    _get_save_callback

Return the Callback that is responsible for saving the expectation values of the system.
=#
function _get_save_callback(sol::AbstractODESolution, method::Type{SF}) where {SF <: AbstractSaveFunc}
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    if hasproperty(kwargs, :callback) && !isnothing(kwargs.callback)
        return _get_save_callback(kwargs.callback, method)
    else
        return nothing
    end
end
_get_save_callback(integrator::AbstractODEIntegrator, method::Type{SF}) where {SF <: AbstractSaveFunc} =
    _get_save_callback(integrator.opts.callback, method)
function _get_save_callback(cb::CallbackSet, method::Type{SF}) where {SF <: AbstractSaveFunc}
    cbs_discrete = cb.discrete_callbacks
    if length(cbs_discrete) > 0
        idx = _get_save_callback_idx(cb, method)
        _cb = cb.discrete_callbacks[idx]
        return _get_save_callback(_cb, method)
    else
        return nothing
    end
end
function _get_save_callback(cb::DiscreteCallback, ::Type{SF}) where {SF <: AbstractSaveFunc}
    if typeof(cb.affect!) <: FunctionCallingAffect && typeof(cb.affect!.func) <: AbstractSaveFunc
        return cb
    end
    return nothing
end
_get_save_callback(cb::ContinuousCallback, ::Type{SF}) where {SF <: AbstractSaveFunc} = nothing

_get_save_callback_idx(cb, method) = 1

# %% ------------ Noise Measurement Helpers ------------ %%

# TODO: To improve. See https://github.com/SciML/DiffEqNoiseProcess.jl/issues/214
function _homodyne_dWdt!(dWdt_cache, integrator, tlist, iter)
    idx = findfirst(>=(tlist[iter[] - 1]), integrator.W.t)

    # We are assuming that the last element is tlist[iter[]]
    @inbounds dWdt_cache .= (integrator.W.u[end] .- integrator.W.u[idx]) ./ (integrator.W.t[end] - integrator.W.t[idx])

    return nothing
end
