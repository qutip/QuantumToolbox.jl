function bloch_redfield_tensor(
        H::QuantumObject{Operator},
        a_ops::Union{AbstractVector, Tuple},
        c_ops::Union{AbstractVector, Tuple, Nothing}=nothing;
        sec_cutoff::Real=0.1,
        fock_basis::Union{Val,Bool}=Val(false)
    )
    rst = eigenstates(H)
    U = QuantumObject(rst.vectors, Operator(), H.dimensions)
    sec_cutoff = float(sec_cutoff)

    # in fock basis
    R0 = liouvillian(H, c_ops)
    
    # set fock_basis=Val(false) and change basis together at the end
    R1 = 0
    if !isempty(a_ops)
        R1 += mapreduce(x -> _brterm(rst, x[1], x[2], sec_cutoff, Val(false)), +, a_ops) 
        
        # do (a_op, spectra)
        #     _brterm(rst, a_op, spectra, sec_cutoff, Val(false))
        # end
    end

    SU = sprepost(U, U') # transformation matrix from eigen basis back to fock basis
    if getVal(fock_basis)
        return R0 + SU * R1 * SU'
    else
        return SU' * R0 * SU + R1, U
    end
end

function brterm(
        H::QuantumObject{Operator},
        a_op::QuantumObject{Operator}, 
        spectra::Function;
        sec_cutoff::Real=0.1, 
        fock_basis::Union{Bool, Val}=Val(false)
    )
    rst = eigenstates(H)
    term = _brterm(rst, a_op, spectra, sec_cutoff, makeVal(fock_basis))
    if getVal(fock_basis)
        return term
    else
        return term, Qobj(rst.vectors, Operator(), rst.dimensions)
    end
end

# method for `EigsolveResult` (so `eigenstates` is not called 
# repeatedly when called by `bloch_redfield_tensor`)
function _brterm(
        rst::EigsolveResult,
        a_op::T,
        spectra::F,
        sec_cutoff::Real,
        fock_basis::Union{Val{true},Val{false}}
    ) where {T<:QuantumObject{Operator},F<:Function}

    _check_br_spectra(spectra)

    U, N = rst.vectors, prod(rst.dimensions) 
    
    skew = @. rst.values - rst.values' |> real
    spectrum = spectra.(skew)
    
    A_mat = U' * a_op.data * U

    ac_term = (A_mat .* spectrum) * A_mat
    bd_term = A_mat * (A_mat .* trans(spectrum))

    if sec_cutoff != -1
        m_cut = similar(skew)
        map!(x -> abs(x) < sec_cutoff, m_cut, skew)
        ac_term .*= m_cut
        bd_term .*= m_cut

        Id = I(N)
        vec_skew = vec(skew)
        M_cut = @. abs(vec_skew - vec_skew') < sec_cutoff
    end
    
    out = 0.5 * (
        + QuantumToolbox._sprepost(A_mat .* trans(spectrum), A_mat)
        + QuantumToolbox._sprepost(A_mat, A_mat .* spectrum)
        - _spost(ac_term, Id) 
        - QuantumToolbox._spre(bd_term, Id)
    )

    (sec_cutoff != -1) && (out .*= M_cut)
    

    if getVal(fock_basis)
        SU = _sprepost(U, U')
        return QuantumObject(SU * out * SU', SuperOperator(), rst.dimensions)
    else
        return QuantumObject(out, SuperOperator(), rst.dimensions)
    end
end

function bloch_redfield_tensor()

end