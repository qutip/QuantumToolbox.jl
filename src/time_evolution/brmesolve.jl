export bloch_redfield_tensor, brterm, brmesolve

@doc raw"""
    bloch_redfield_tensor(
        H::QuantumObject{Operator},
        a_ops::Union{AbstractVector, Tuple, Nothing},
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
        a_ops::Union{AbstractVector, Tuple, Nothing},
        c_ops::Union{AbstractVector, Tuple, Nothing} = nothing;
        sec_cutoff::Real = 0.1,
        fock_basis::Union{Val, Bool} = Val(false),
    )
    rst = eigenstates(H)
    U = QuantumObject(rst.vectors, Operator(), H.dimensions)
    sec_cutoff = float(sec_cutoff)

    H_new = getVal(fock_basis) ? H : QuantumObject(Diagonal(rst.values), Operator(), H.dimensions)
    c_ops_new = isnothing(c_ops) ? nothing : (getVal(fock_basis) ? c_ops : map(x -> to_sparse_if_needed(Val(issparse(x)), U' * x * U), c_ops))
    L0 = liouvillian(H_new, c_ops_new)

    # Check whether we can rotate the terms to the eigenbasis directly in the Hamiltonian space
    fock_basis_hamiltonian = getVal(fock_basis) && sec_cutoff == -1

    R =
        (isnothing(a_ops) || isempty(a_ops)) ? 0 :
        sum(x -> _brterm(rst, x[1], x[2], sec_cutoff, fock_basis_hamiltonian), a_ops)

    # If in fock basis, we need to transform the terms back to the fock basis
    # Note: we can transform the terms in the Hamiltonian space only if sec_cutoff is -1
    # otherwise, we need to use the SU superoperator below to transform the entire Liouvillian
    # at the end, due to the action of M_cut
    if getVal(fock_basis)
        if fock_basis_hamiltonian
            return L0 + R # Already rotated in the Hamiltonian space
        else
            SU = sprepost(U, U')
            return L0 + SU * R * SU'
        end
    else
        return L0 + R, U
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
        fock_basis::Union{Bool, Val} = Val(false),
    )
    rst = eigenstates(H)
    U = QuantumObject(rst.vectors, Operator(), H.dimensions)

    # Check whether we can rotate the terms to the eigenbasis directly in the Hamiltonian space
    fock_basis_hamiltonian = getVal(fock_basis) && sec_cutoff == -1

    term = _brterm(rst, a_op, spectra, sec_cutoff, fock_basis_hamiltonian)
    if getVal(fock_basis)
        if fock_basis_hamiltonian
            return term # Already rotated in the Hamiltonian space
        else
            SU = sprepost(U, U')
            return SU * term * SU'
        end
    else
        return term, U
    end
end

function _brterm(
        rst::EigsolveResult,
        a_op::T,
        spectra::F,
        sec_cutoff::Real,
        fock_basis_hamiltonian::Union{Bool, Val},
    ) where {T <: QuantumObject{Operator}, F <: Function}
    _check_br_spectra(spectra)

    U = rst.vectors

    skew = @. rst.values - rst.values' |> real
    spectrum = spectra.(skew)

    A_mat = to_sparse_if_needed(Val(issparse(a_op.data)), U' * a_op.data * U)
    A_mat_spec = A_mat .* spectrum
    A_mat_spec_t = A_mat .* transpose(spectrum)

    ac_term = A_mat_spec * A_mat
    bd_term = A_mat * A_mat_spec_t

    if sec_cutoff != -1
        m_cut = map(x -> abs(x) < sec_cutoff, skew)
        ac_term .*= m_cut
        bd_term .*= m_cut
    end

    # Rotate the terms to the eigenbasis if possible
    if getVal(fock_basis_hamiltonian)
        A_mat = to_sparse_if_needed(Val(issparse(A_mat)), U * A_mat * U')
        A_mat_spec = to_sparse_if_needed(Val(issparse(A_mat_spec)), U * A_mat_spec * U')
        A_mat_spec_t = to_sparse_if_needed(Val(issparse(A_mat_spec_t)), U * A_mat_spec_t * U')
        ac_term = to_sparse_if_needed(Val(issparse(ac_term)), U * ac_term * U')
        bd_term = to_sparse_if_needed(Val(issparse(bd_term)), U * bd_term * U')
    end

    # Remove small values before passing in the Liouville space
    if settings.auto_tidyup
        tidyup!(A_mat)
        tidyup!(A_mat_spec)
        tidyup!(A_mat_spec_t)
        tidyup!(ac_term)
        tidyup!(bd_term)
    end

    out = (_sprepost(A_mat_spec_t, A_mat) + _sprepost(A_mat, A_mat_spec) - _spost(ac_term) - _spre(bd_term)) / 2

    if (sec_cutoff != -1)
        vec_skew = vec(skew)
        M_cut = @. abs(vec_skew - vec_skew') < sec_cutoff

        out .*= M_cut
    end

    return QuantumObject(out, SuperOperator(), rst.dimensions)
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
        a_ops::Union{Nothing, AbstractVector, Tuple},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        sec_cutoff::Real = 0.1,
        e_ops::Union{Nothing, AbstractVector} = nothing,
        kwargs...,
    )
    R = bloch_redfield_tensor(H, a_ops, c_ops; sec_cutoff = sec_cutoff, fock_basis = Val(true))

    return mesolve(R, ψ0, tlist, nothing; e_ops = e_ops, kwargs...)
end

function _check_br_spectra(f::Function)
    length(methods(f, [Real]).ms) == 0 &&
        throw(ArgumentError("The following function must only accept one argument: `$(nameof(f))(ω)` with ω<:Real"))
    return nothing
end
