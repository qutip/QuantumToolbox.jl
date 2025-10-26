#=
Helper functions for the mesolve callbacks.
=#

struct SaveFuncMESolve{TE,PT<:Union{Nothing,Progress},IT,TEXPV<:Union{Nothing,AbstractMatrix}} <: AbstractSaveFunc
    e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncMESolve)(u, t, integrator) = _save_func_mesolve(u, integrator, f.e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncMESolve{Nothing})(u, t, integrator) = _save_func(integrator, f.progr)

_get_e_ops_data(e_ops, ::Type{SaveFuncMESolve}) = [_generate_mesolve_e_op(op) for op in e_ops] # Broadcasting generates type instabilities on Julia v1.10

_get_progress_desc(::Type{SaveFuncMESolve}) = "[mesolve] "

##

# When e_ops is a list of operators
function _save_func_mesolve(u, integrator, e_ops, progr, iter, expvals)
    # This is equivalent to tr(op * ρ), when both are matrices.
    # The advantage of using this convention is that We don't need
    # to reshape u to make it a matrix, but we reshape the e_ops once.

    ρ = u
    _expect = op -> dot(op, ρ)
    @. expvals[:, iter[]] = _expect(e_ops)
    iter[] += 1

    return _save_func(integrator, progr)
end

function _mesolve_callbacks_new_e_ops!(integrator::AbstractODEIntegrator, e_ops)
    cb = _get_save_callback(integrator, SaveFuncMESolve)
    if cb isa Nothing
        return nothing
    else
        cb.affect!.func.e_ops .= e_ops # Only works if e_ops is a Vector of operators
        return nothing
    end
end

_generate_mesolve_e_op(op) = mat2vec(adjoint(get_data(op)))
