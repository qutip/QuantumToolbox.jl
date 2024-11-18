#=
Helper functions for the sesolve callbacks.
=#

struct SaveFuncSESolve{TE,PT<:Union{Nothing,ProgressBar},IT,TEXPV<:Union{Nothing,AbstractMatrix}}
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncSESolve)(integrator) = _save_func_sesolve(integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncSESolve{Nothing})(integrator) = _save_func(integrator, f.progr) # Common for both mesolve and sesolve

##

# When e_ops is a list of operators
function _save_func_sesolve(integrator, e_ops, progr, iter, expvals)
    ψ = integrator.u
    _expect = op -> dot(ψ, op, ψ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end
