#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

##

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

_get_e_ops_data(e_ops, ::Type{SaveFuncSESolve}) = get_data.(e_ops)
_get_e_ops_data(e_ops, ::Type{SaveFuncMESolve}) = [_generate_mesolve_e_op(op) for op in e_ops] # Broadcasting generates type instabilities on Julia v1.10

_generate_mesolve_e_op(op) = mat2vec(adjoint(get_data(op)))

##

# When e_ops is Nothing. Common for both mesolve and sesolve
function _save_func(integrator, progr)
    next!(progr)
    u_modified!(integrator, false)
    return nothing
end

# When progr is Nothing. Common for both mesolve and sesolve
function _save_func(integrator, progr::Nothing)
    u_modified!(integrator, false)
    return nothing
end

##

# Get the e_ops from a given AbstractODESolution. Valid for `sesolve`, `mesolve` and `ssesolve`.
function _se_me_sse_get_expvals(sol::AbstractODESolution)
    cb = _se_me_sse_get_save_callback(sol)
    if cb isa Nothing
        return nothing
    else
        return cb.affect!.expvals
    end
end

function _se_me_sse_get_save_callback(sol::AbstractODESolution)
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    if hasproperty(kwargs, :callback)
        return _se_me_sse_get_save_callback(kwargs.callback)
    else
        return nothing
    end
end
_se_me_sse_get_save_callback(integrator::AbstractODEIntegrator) = _se_me_sse_get_save_callback(integrator.opts.callback)
function _se_me_sse_get_save_callback(cb::CallbackSet)
    cbs_discrete = cb.discrete_callbacks
    if length(cbs_discrete) > 0
        _cb = cb.discrete_callbacks[1]
        return _se_me_sse_get_save_callback(_cb)
    else
        return nothing
    end
end
_se_me_sse_get_save_callback(cb::DiscreteCallback) =
    if (cb.affect! isa SaveFuncSESolve) || (cb.affect! isa SaveFuncMESolve)
        return cb
    else
        return nothing
    end
_se_me_sse_get_save_callback(cb::ContinuousCallback) = nothing
