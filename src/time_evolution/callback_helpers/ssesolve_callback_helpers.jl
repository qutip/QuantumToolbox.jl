#=
Helper functions for the ssesolve callbacks. Almost equal to the sesolve case, but with an additional possibility to store the measurement operators expectation values. Also, this callback is not the first one, but the second one, due to the presence of the normalization callback.
=#

struct SaveFuncSSESolve{
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

(f::SaveFuncSSESolve)(u, t, integrator) =
    _save_func_ssesolve(u, integrator, f.e_ops, f.m_ops, f.progr, f.iter, f.expvals, f.m_expvals)
(f::SaveFuncSSESolve{false,Nothing})(u, t, integrator) = _save_func(integrator, f.progr) # Common for both all solvers

_get_e_ops_data(e_ops, ::Type{SaveFuncSSESolve}) = get_data.(e_ops)
_get_m_ops_data(sc_ops, ::Type{SaveFuncSSESolve}) = map(op -> Hermitian(get_data(op) + get_data(op)'), sc_ops)

_get_save_callback_idx(cb, ::Type{SaveFuncSSESolve}) = 2 # The first one is the normalization callback

##

# When e_ops is a list of operators
function _save_func_ssesolve(u, integrator, e_ops, m_ops, progr, iter, expvals, m_expvals)
    ψ = u

    _expect = op -> dot(ψ, op, ψ)

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

##

#=
    This function generates the normalization callback. It is needed to stabilize the integration when using the ssesolve method.
=#
function _ssesolve_generate_normalize_cb()
    _condition = (u, t, integrator) -> true
    _affect! = (integrator) -> normalize!(integrator.u)
    cb = DiscreteCallback(_condition, _affect!; save_positions = (false, false))

    return cb
end
