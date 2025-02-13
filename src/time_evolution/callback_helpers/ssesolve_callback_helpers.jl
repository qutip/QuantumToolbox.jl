#=
Helper functions for the ssesolve callbacks. Equal to the sesolve case, but with an additional normalization before saving the expectation values.
=#

struct SaveFuncSSESolve{TE,PT<:Union{Nothing,ProgressBar},IT,TEXPV<:Union{Nothing,AbstractMatrix}} <: AbstractSaveFunc
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncSSESolve)(integrator) = _save_func_ssesolve(integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncSSESolve{Nothing})(integrator) = _save_func(integrator, f.progr) # Common for both all solvers

_get_e_ops_data(e_ops, ::Type{SaveFuncSSESolve}) = get_data.(e_ops)

_get_save_callback_idx(cb, ::Type{SaveFuncSSESolve}) = 1 # The first one is the normalization callback

##

# When e_ops is a list of operators
function _save_func_ssesolve(integrator, e_ops, progr, iter, expvals)
    ψ = normalize!(integrator.u)
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end

##

#=
This function adds the normalization callback to the kwargs. It is needed to stabilize the integration when using the ssesolve method.
=#
function _ssesolve_add_normalize_cb(kwargs)
    _condition = (u, t, integrator) -> true
    _affect! = (integrator) -> normalize!(integrator.u)
    cb = DiscreteCallback(_condition, _affect!; save_positions = (false, false))

    cb_set = haskey(kwargs, :callback) ? CallbackSet(kwargs[:callback], cb) : cb

    kwargs2 = merge(kwargs, (callback = cb_set,))

    return kwargs2
end
