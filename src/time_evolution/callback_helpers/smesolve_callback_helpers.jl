#=
Helper functions for the smesolve callbacks. Almost equal to the mesolve case, but with an additional possibility to store the measurement operators expectation values. 
=#

struct SaveFuncSMESolve{
    SM,
    TE,
    TE2,
    PT<:Union{Nothing,ProgressBar},
    IT,
    TEXPV<:Union{Nothing,AbstractMatrix},
    TMEXPV<:Union{Nothing,AbstractMatrix},
} <: AbstractSaveFunc
    store_measurement::Val{SM}
    e_ops::TE
    m_ops::TE2
    progr::PT
    iter::IT
    expvals::TEXPV
    m_expvals::TMEXPV
end

(f::SaveFuncSMESolve)(integrator) =
    _save_func_smesolve(integrator, f.e_ops, f.m_ops, f.progr, f.iter, f.expvals, f.m_expvals)
(f::SaveFuncSMESolve{false,Nothing})(integrator) = _save_func(integrator, f.progr) # Common for both all solvers

_get_e_ops_data(e_ops, ::Type{SaveFuncSMESolve}) = _get_e_ops_data(e_ops, SaveFuncMESolve)
_get_m_ops_data(sc_ops, ::Type{SaveFuncSMESolve}) =
    map(op -> _generate_mesolve_e_op(op) + _generate_mesolve_e_op(op'), sc_ops)

##

# When e_ops is a list of operators
function _save_func_smesolve(integrator, e_ops, m_ops, progr, iter, expvals, m_expvals)
    # This is equivalent to tr(op * ρ), when both are matrices.
    # The advantage of using this convention is that We don't need
    # to reshape u to make it a matrix, but we reshape the e_ops once.

    ρ = integrator.u

    _expect = op -> dot(op, ρ)

    if !isnothing(e_ops)
        @. expvals[:, iter[]] = _expect(e_ops)
    end

    if !isnothing(m_expvals) && iter[] > 1
        _dWdt = _homodyne_dWdt(integrator)
        @. m_expvals[:, iter[]-1] = real(_expect(m_ops)) + _dWdt
    end

    iter[] += 1

    return _save_func(integrator, progr)
end
