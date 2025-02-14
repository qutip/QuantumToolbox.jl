#=
Helper functions for the ssesolve callbacks. Equal to the sesolve case, but with an additional normalization before saving the expectation values.
=#

struct SaveFuncSSESolve{
    SM,
    TE,
    TME,
    PT<:Union{Nothing,ProgressBar},
    IT,
    TEXPV<:Union{Nothing,AbstractMatrix},
    TMEXPV<:Union{Nothing,AbstractMatrix},
} <: AbstractSaveFunc
    store_measurement::Val{SM}
    e_ops::TE
    m_ops::TME
    progr::PT
    iter::IT
    expvals::TEXPV
    m_expvals::TMEXPV
end

(f::SaveFuncSSESolve)(integrator) =
    _save_func_ssesolve(integrator, f.e_ops, f.m_ops, f.progr, f.iter, f.expvals, f.m_expvals)
(f::SaveFuncSSESolve{false,Nothing})(integrator) = _save_func(integrator, f.progr) # Common for both all solvers

_get_e_ops_data(e_ops, ::Type{SaveFuncSSESolve}) = get_data.(e_ops)
_get_m_ops_data(sc_ops, ::Type{SaveFuncSSESolve}) = map(op -> Hermitian(get_data(op) + get_data(op)'), sc_ops)

_get_save_callback_idx(cb, ::Type{SaveFuncSSESolve}) = 2 # The first one is the normalization callback

##

# When e_ops is a list of operators
function _save_func_ssesolve(integrator, e_ops, m_ops, progr, iter, expvals, m_expvals)
    ψ = integrator.u

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

# TODO: Add some cache mechanism to avoid memory allocations
function _homodyne_dWdt(integrator)
    @inbounds _dWdt = (integrator.W.u[end] .- integrator.W.u[end-1]) ./ (integrator.W.t[end] - integrator.W.t[end-1])

    return _dWdt
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
