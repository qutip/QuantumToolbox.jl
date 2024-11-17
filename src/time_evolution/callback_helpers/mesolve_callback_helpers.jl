#=
Helper functions for the mesolve callbacks.
=#

struct SaveFuncMESolve{TE,PT<:Union{Nothing,ProgressBar},IT,TEXPV<:Union{Nothing,AbstractMatrix}}
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncMESolve)(integrator) = _save_func_mesolve(integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncMESolve{Nothing})(integrator) = _save_func(integrator, f.progr)

##

# When e_ops is a list of operators
function _save_func_mesolve(integrator, e_ops, progr, iter, expvals)
    # This is equivalent to tr(op * ρ), when both are matrices.
    # The advantage of using this convention is that We don't need
    # to reshape u to make it a matrix, but we reshape the e_ops once.

    ρ = integrator.u
    _expect = op -> dot(op, ρ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end
