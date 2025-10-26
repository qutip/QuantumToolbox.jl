#=
Helper functions for the sesolve callbacks.
=#

struct SaveFuncSESolve{TE,PT<:Union{Nothing,Progress},IT,TEXPV<:Union{Nothing,AbstractMatrix}} <: AbstractSaveFunc
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncSESolve)(u, t, integrator) = _save_func_sesolve(u, integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncSESolve{Nothing})(u, t, integrator) = _save_func(integrator, f.progr) # Common for both mesolve and sesolve

_get_e_ops_data(e_ops, ::Type{SaveFuncSESolve}) = get_data.(e_ops)

_get_progress_desc(::Type{SaveFuncSESolve}) = "[sesolve] "

##

# When e_ops is a list of operators
function _save_func_sesolve(u, integrator, e_ops, progr, iter, expvals)
    ψ = u
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end
