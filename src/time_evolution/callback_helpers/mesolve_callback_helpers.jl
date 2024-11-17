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

function _mesolve_callbacks_new_e_ops!(integrator::AbstractODEIntegrator, e_ops)
    cb = _se_me_sse_get_save_callback(integrator)
    if cb isa Nothing
        return nothing
    else
        cb.affect!.e_ops .= e_ops # Only works if e_ops is a Vector of operators
        return nothing
    end
end
