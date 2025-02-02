#=
Helper functions for the ssesolve callbacks. Equal to the sesolve case, but with an additional normalization before saving the expectation values.
=#

struct SaveFuncSSESolve{TE,PT<:Union{Nothing,ProgressBar},IT,TEXPV<:Union{Nothing,AbstractMatrix}}
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncSSESolve)(integrator) = _save_func_ssesolve(integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncSSESolve{Nothing})(integrator) = _save_func(integrator, f.progr) # Common for both mesolve and sesolve

##

# When e_ops is a list of operators
function _save_func_ssesolve(integrator, e_ops, progr, iter, expvals)
    ψ = normalize!(integrator.u)
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end
