export bloch_redfield_tensor, brterm, brmesolve

@doc raw"""
    bloch_redfield_tensor(
        H::QuantumObject{Operator},
        a_ops::Union{AbstractVector, Tuple},
        c_ops::Union{AbstractVector, Tuple, Nothing}=nothing;
        sec_cutoff::Real=0.1,
        fock_basis::Union{Val,Bool}=Val(false)
    )

Calculates the Bloch-Redfield tensor ([`SuperOperator`](@ref)) for a system given a set of operators and corresponding spectral functions that describes the system's coupling to its environment.

## Arguments

- `H`: The system Hamiltonian. Must be an [`Operator`](@ref)
- `a_ops`: Nested list with each element is a `Tuple` of operator-function pairs `(a_op, spectra)`, and the coupling [`Operator`](@ref) `a_op` must be hermitian with corresponding `spectra` being a `Function` of transition energy
- `c_ops`: List of collapse operators corresponding to Lindblad dissipators
- `sec_cutoff`: Cutoff for secular approximation. Use `-1` if secular approximation is not used when evaluating bath-coupling terms.
- `fock_basis`: Whether to return the tensor in the input (fock) basis or the diagonalized (eigen) basis.

## Return

The return depends on `fock_basis`.

- `fock_basis=Val(true)`: return the Bloch-Redfield tensor (in the fock basis) only.
- `fock_basis=Val(false)`: return the Bloch-Redfield tensor (in the eigen basis) along with the transformation matrix from eigen to fock basis.
"""
function bloch_redfield_tensor(
    H::QuantumObject{Operator},
    a_ops::Union{AbstractVector,Tuple},
    c_ops::Union{AbstractVector,Tuple,Nothing} = nothing;
    sec_cutoff::Real = 0.1,
    fock_basis::Union{Val,Bool} = Val(false),
)
    rst = eigenstates(H)
    U = QuantumObject(rst.vectors, Operator(), H.dimensions)
    sec_cutoff = float(sec_cutoff)

    # in fock basis
    R0 = liouvillian(H, c_ops)

    # set fock_basis=Val(false) and change basis together at the end
    R1 = 0
    isempty(a_ops) || (R1 += mapreduce(x -> _brterm(rst, x[1], x[2], sec_cutoff, Val(false)), +, a_ops))

    SU = sprepost(U, U') # transformation matrix from eigen basis back to fock basis
    if getVal(fock_basis)
        return R0 + SU * R1 * SU'
    else
        return SU' * R0 * SU + R1, U
    end
end

@doc raw"""
    brterm(
        H::QuantumObject{Operator},
        a_op::QuantumObject{Operator}, 
        spectra::Function;
        sec_cutoff::Real=0.1, 
        fock_basis::Union{Bool, Val}=Val(false)
    )

Calculates the contribution of one coupling operator to the Bloch-Redfield tensor.

## Argument

- `H`: The system Hamiltonian. Must be an [`Operator`](@ref)
- `a_op`: The operator coupling to the environment. Must be hermitian.
- `spectra`: The corresponding environment spectra as a `Function` of transition energy.
- `sec_cutoff`: Cutoff for secular approximation. Use `-1` if secular approximation is not used when evaluating bath-coupling terms.
- `fock_basis`: Whether to return the tensor in the input (fock) basis or the diagonalized (eigen) basis.

## Return

The return depends on `fock_basis`.

- `fock_basis=Val(true)`: return the Bloch-Redfield term (in the fock basis) only.
- `fock_basis=Val(false)`: return the Bloch-Redfield term (in the eigen basis) along with the transformation matrix from eigen to fock basis.
"""
function brterm(
    H::QuantumObject{Operator},
    a_op::QuantumObject{Operator},
    spectra::Function;
    sec_cutoff::Real = 0.1,
    fock_basis::Union{Bool,Val} = Val(false),
)
    rst = eigenstates(H)
    term = _brterm(rst, a_op, spectra, sec_cutoff, makeVal(fock_basis))
    if getVal(fock_basis)
        return term
    else
        return term, Qobj(rst.vectors, Operator(), rst.dimensions)
    end
end

function _brterm(
    rst::EigsolveResult,
    a_op::T,
    spectra::F,
    sec_cutoff::Real,
    fock_basis::Union{Val{true},Val{false}},
) where {T<:QuantumObject{Operator},F<:Function}
    _check_br_spectra(spectra)

    U = rst.vectors
    Id = I(prod(rst.dimensions))

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

        vec_skew = vec(skew)
        M_cut = @. abs(vec_skew - vec_skew') < sec_cutoff
    end

    out =
        0.5 * (
            + _sprepost(A_mat .* trans(spectrum), A_mat) + _sprepost(A_mat, A_mat .* spectrum) - _spost(ac_term, Id) -
            _spre(bd_term, Id)
        )

    (sec_cutoff != -1) && (out .*= M_cut)

    if getVal(fock_basis)
        SU = _sprepost(U, U')
        return QuantumObject(SU * out * SU', SuperOperator(), rst.dimensions)
    else
        return QuantumObject(out, SuperOperator(), rst.dimensions)
    end
end

@doc raw"""
    brmesolve(
        H::QuantumObject{Operator}, 
        ψ0::QuantumObject,
        tlist::AbstractVector, 
        a_ops::Union{Nothing, AbstractVector, Tuple}, 
        c_ops::Union{Nothing, AbstractVector, Tuple}=nothing;
        sec_cutoff::Real=0.1,
        e_ops::Union{Nothing, AbstractVector}=nothing,
        kwargs...,
    )

Solves for the dynamics of a system using the Bloch-Redfield master equation, given an input Hamiltonian, Hermitian bath-coupling terms and their associated spectral functions, as well as possible Lindblad collapse operators.

## Arguments

- `H`: The system Hamiltonian. Must be an [`Operator`](@ref)
- `ψ0`:  Initial state of the system $|\psi(0)\rangle$. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `a_ops`: Nested list with each element is a `Tuple` of operator-function pairs `(a_op, spectra)`, and the coupling [`Operator`](@ref) `a_op` must be hermitian with corresponding `spectra` being a `Function` of transition energy
- `c_ops`: List of collapse operators corresponding to Lindblad dissipators
- `sec_cutoff`: Cutoff for secular approximation. Use `-1` if secular approximation is not used when evaluating bath-coupling terms.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `kwargs`: Keyword arguments for [`mesolve`](@ref).

## Notes

- This function will automatically generate the [`bloch_redfield_tensor`](@ref) and solve the time evolution with [`mesolve`](@ref).

# Returns

- `sol::TimeEvolutionSol`: The solution of the time evolution. See also [`TimeEvolutionSol`](@ref)
"""
function brmesolve(
    H::QuantumObject{Operator},
    ψ0::QuantumObject,
    tlist::AbstractVector,
    a_ops::Union{Nothing,AbstractVector,Tuple},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    sec_cutoff::Real = 0.1,
    e_ops::Union{Nothing,AbstractVector} = nothing,
    kwargs...,
)
    R = bloch_redfield_tensor(H, a_ops, c_ops; sec_cutoff = sec_cutoff, fock_basis = Val(true))

    return mesolve(R, ψ0, tlist, nothing; e_ops = e_ops, kwargs...)
end

function _check_br_spectra(f::Function)
    meths = methods(f, [Real])
    length(meths.ms) == 0 &&
        throw(ArgumentError("The following function must accept one argument: `$(meths.mt.name)(ω)` with ω<:Real"))
    return nothing
end
