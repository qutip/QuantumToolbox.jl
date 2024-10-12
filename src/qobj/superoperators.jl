#=
Functions for generating (common) quantum super-operators.
=#

export spre, spost, sprepost, lindblad_dissipator

# intrinsic functions for super-operators
# (keep these because they take AbstractMatrix as input)
_spre(A::AbstractMatrix, Id::AbstractMatrix) = kron(Id, sparse(A))
_spre(A::AbstractSparseMatrix, Id::AbstractMatrix) = kron(Id, A)
_spost(B::AbstractMatrix, Id::AbstractMatrix) = kron(transpose(sparse(B)), Id)
_spost(B::AbstractSparseMatrix, Id::AbstractMatrix) = kron(transpose(B), Id)
_sprepost(A::AbstractMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), sparse(A))
_sprepost(A::AbstractMatrix, B::AbstractSparseMatrix) = kron(transpose(B), sparse(A))
_sprepost(A::AbstractSparseMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), A)
_sprepost(A::AbstractSparseMatrix, B::AbstractSparseMatrix) = kron(transpose(B), A)

@doc raw"""
    spre(A::QuantumObject, Id_cache=I(size(A,1)))

Returns the [`SuperOperator`](@ref) form of `A` acting on the left of the density matrix operator: ``\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{\mathbb{1}} \otimes \hat{A}``, namely 

```math
\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{\mathbb{1}} \otimes \hat{A} ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
spre(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, Id_cache = I(size(A, 1))) where {T} =
    QuantumObject(_spre(A.data, Id_cache), SuperOperator, A.dims)

@doc raw"""
    spost(B::QuantumObject, Id_cache=I(size(B,1)))

Returns the [`SuperOperator`](@ref) form of `B` acting on the right of the density matrix operator: ``\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{\mathbb{1}}``, namely

```math
\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{\mathbb{1}} ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
spost(B::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, Id_cache = I(size(B, 1))) where {T} =
    QuantumObject(_spost(B.data, Id_cache), SuperOperator, B.dims)

@doc raw"""
    sprepost(A::QuantumObject, B::QuantumObject)

Returns the [`SuperOperator`](@ref) form of `A` and `B` acting on the left and right of the density matrix operator, respectively: ``\mathcal{O} \left( \hat{A}, \hat{B} \right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{A}``, namely

```math
\mathcal{O} \left(\hat{A}, \hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{A} ~ |\hat{\rho}\rangle\rangle = \textrm{spre}(\hat{A}) * \textrm{spost}(\hat{B}) ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

See also [`spre`](@ref) and [`spost`](@ref).
"""
function sprepost(
    A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects don't have the same Hilbert dimension."))

    return QuantumObject(_sprepost(A.data, B.data), SuperOperator, A.dims)
end

@doc raw"""
    lindblad_dissipator(O::QuantumObject, Id_cache=I(size(O,1))

Returns the Lindblad [`SuperOperator`](@ref) defined as

```math
\mathcal{D} \left( \hat{O} \right) \left[ \hat{\rho} \right] = \frac{1}{2} \left( 2 \hat{O} \hat{\rho} \hat{O}^\dagger - 
\hat{O}^\dagger \hat{O} \hat{\rho} - \hat{\rho} \hat{O}^\dagger \hat{O} \right)
```

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.

See also [`spre`](@ref) and [`spost`](@ref).
"""
function lindblad_dissipator(
    O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject},
    Id_cache = I(size(O, 1)),
) where {T}
    Od_O = O' * O
    return sprepost(O, O') - spre(Od_O, Id_cache) / 2 - spost(Od_O, Id_cache) / 2
end

# It is already a SuperOperator
lindblad_dissipator(O::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject}, Id_cache) where {T} = O
