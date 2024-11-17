#=
This file contains helper functions for callbacks. The affect! function are defined taking advantage of the Julia struct, which allows to store some cache exclusively for the callback.
=#

include("sesolve_callback_helpers.jl")
include("mesolve_callback_helpers.jl")
include("mcsolve_callback_helpers.jl")

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
    e_ops_data = e_ops isa Nothing ? nothing : get_data.(e_ops)

    progr = getVal(progress_bar) ? ProgressBar(length(tlist), enable = getVal(progress_bar)) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))

    _save_affect! = method(e_ops_data, progr, Ref(1), expvals)
    return PresetTimeCallback(tlist, _save_affect!, save_positions = (false, false))
end

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
    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple to support Zygote.jl
    if hasproperty(kwargs, :callback)
        return _se_me_sse_get_expvals(kwargs.callback)
    else
        return nothing
    end
end
function _se_me_sse_get_expvals(cb::CallbackSet)
    _cb = cb.discrete_callbacks[1]
    return _se_me_sse_get_expvals(_cb)
end
_se_me_sse_get_expvals(cb::DiscreteCallback) =
    if (cb.affect! isa SaveFuncSESolve) || (cb.affect! isa SaveFuncMESolve)
        return cb.affect!.expvals
    else
        return nothing
    end
_se_me_sse_get_expvals(cb::ContinuousCallback) = nothing
