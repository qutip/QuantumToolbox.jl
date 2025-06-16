#=
Functions for generating (common) quantum super-operators.
=#

export spre, spost, sprepost, liouvillian, lindblad_dissipator

# intrinsic functions for super-operators
## keep these because they take AbstractMatrix as input and ensure the output is sparse matrix
_spre(A::AbstractMatrix, Id::AbstractMatrix) = kron(Id, sparse(A))
_spre(A::AbstractSparseMatrix, Id::AbstractMatrix) = kron(Id, A)
_spost(B::AbstractMatrix, Id::AbstractMatrix) = kron(transpose(sparse(B)), Id)
_spost(B::AbstractSparseMatrix, Id::AbstractMatrix) = kron(transpose(B), Id)
_sprepost(A::AbstractMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), sparse(A))
_sprepost(A::AbstractMatrix, B::AbstractSparseMatrix) = kron(transpose(B), sparse(A))
_sprepost(A::AbstractSparseMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), A)
_sprepost(A::AbstractSparseMatrix, B::AbstractSparseMatrix) = kron(transpose(B), A)
function _sprepost(A, B) # for any other input types
    # TODO: use the commented code (since it is optimized for certain types of SciMLOperators, and was able to give correct results before `SciMLOperators v1.0.0`)
    # Id_cache = I(size(A, 1))
    # return _spre(A, Id_cache) * _spost(B, Id_cache)

    return kron(transpose(B), A)
end

## if input is AbstractSciMLOperator 
## some of them are optimized to speed things up
## the rest of the SciMLOperators will just use lazy tensor (and prompt a warning)
_spre(A::MatrixOperator, Id::AbstractMatrix) = MatrixOperator(_spre(A.A, Id))
_spre(A::ScaledOperator, Id::AbstractMatrix) = ScaledOperator(A.λ, _spre(A.L, Id))
_spre(A::AddedOperator, Id::AbstractMatrix) = AddedOperator(map(op -> _spre(op, Id), A.ops))
function _spre(A::AbstractSciMLOperator, Id::AbstractMatrix)
    _lazy_tensor_warning(Id, A)
    return kron(Id, A)
end

_spost(B::MatrixOperator, Id::AbstractMatrix) = MatrixOperator(_spost(B.A, Id))
_spost(B::ScaledOperator, Id::AbstractMatrix) = ScaledOperator(B.λ, _spost(B.L, Id))
_spost(B::AddedOperator, Id::AbstractMatrix) = AddedOperator(map(op -> _spost(op, Id), B.ops))
function _spost(B::AbstractSciMLOperator, Id::AbstractMatrix)
    B_T = transpose(B)
    _lazy_tensor_warning(B_T, Id)
    return kron(B_T, Id)
end

## intrinsic liouvillian 
_liouvillian(H::MT, Id::AbstractMatrix) where {MT<:Union{AbstractMatrix,AbstractSciMLOperator}} =
    -1im * (_spre(H, Id) - _spost(H, Id))
_liouvillian(H::MatrixOperator, Id::AbstractMatrix) = MatrixOperator(_liouvillian(H.A, Id))
_liouvillian(H::ScaledOperator, Id::AbstractMatrix) = ScaledOperator(H.λ, _liouvillian(H.L, Id))
_liouvillian(H::AddedOperator, Id::AbstractMatrix) = AddedOperator(map(op -> _liouvillian(op, Id), H.ops))

# intrinsic lindblad_dissipator
function _lindblad_dissipator(O::MT, Id::AbstractMatrix) where {MT<:Union{AbstractMatrix,AbstractSciMLOperator}}
    Od_O = O' * O
    return _sprepost(O, O') - (_spre(Od_O, Id) + _spost(Od_O, Id)) / 2
end
function _lindblad_dissipator(O::MatrixOperator, Id::AbstractMatrix)
    _O = O.A
    Od_O = _O' * _O
    return MatrixOperator(_sprepost(_O, _O') - (_spre(Od_O, Id) + _spost(Od_O, Id)) / 2)
end
function _lindblad_dissipator(O::ScaledOperator, Id::AbstractMatrix)
    λc_λ = conj(O.λ) * O.λ
    return ScaledOperator(λc_λ, _lindblad_dissipator(O.L, Id))
end

@doc raw"""
    spre(A::AbstractQuantumObject, Id_cache=I(size(A,1)))

Returns the [`SuperOperator`](@ref) form of `A` acting on the left of the density matrix operator: ``\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{\mathbb{1}} \otimes \hat{A}``, namely 

```math
\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{\mathbb{1}} \otimes \hat{A} ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when the same function is applied multiple times with a known Hilbert space dimension.

See also [`spost`](@ref) and [`sprepost`](@ref).
"""
spre(A::AbstractQuantumObject{Operator}, Id_cache = I(size(A, 1))) =
    get_typename_wrapper(A)(_spre(A.data, Id_cache), SuperOperator(), A.dimensions)

@doc raw"""
    spost(B::AbstractQuantumObject, Id_cache=I(size(B,1)))

Returns the [`SuperOperator`](@ref) form of `B` acting on the right of the density matrix operator: ``\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{\mathbb{1}}``, namely

```math
\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{\mathbb{1}} ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when the same function is applied multiple times with a known Hilbert space dimension.

See also [`spre`](@ref) and [`sprepost`](@ref).
"""
spost(B::AbstractQuantumObject{Operator}, Id_cache = I(size(B, 1))) =
    get_typename_wrapper(B)(_spost(B.data, Id_cache), SuperOperator(), B.dimensions)

@doc raw"""
    sprepost(A::AbstractQuantumObject, B::AbstractQuantumObject)

Returns the [`SuperOperator`](@ref) form of `A` and `B` acting on the left and right of the density matrix operator, respectively: ``\mathcal{O} \left( \hat{A}, \hat{B} \right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{A}``, namely

```math
\mathcal{O} \left(\hat{A}, \hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{A} ~ |\hat{\rho}\rangle\rangle = \textrm{spre}(\hat{A}) * \textrm{spost}(\hat{B}) ~ |\hat{\rho}\rangle\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

See also [`spre`](@ref) and [`spost`](@ref).
"""
function sprepost(A::AbstractQuantumObject{Operator}, B::AbstractQuantumObject{Operator})
    check_dimensions(A, B)
    return promote_op_type(A, B)(_sprepost(A.data, B.data), SuperOperator(), A.dimensions)
end

@doc raw"""
    lindblad_dissipator(O::AbstractQuantumObject, Id_cache=I(size(O,1))

Returns the Lindblad [`SuperOperator`](@ref) defined as

```math
\mathcal{D} \left( \hat{O} \right) \left[ \hat{\rho} \right] = \frac{1}{2} \left( 2 \hat{O} \hat{\rho} \hat{O}^\dagger - 
\hat{O}^\dagger \hat{O} \hat{\rho} - \hat{\rho} \hat{O}^\dagger \hat{O} \right)
```

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when the same function is applied multiple times with a known Hilbert space dimension.

See also [`spre`](@ref), [`spost`](@ref), and [`sprepost`](@ref).
"""
lindblad_dissipator(O::AbstractQuantumObject{Operator}, Id_cache = I(size(O, 1))) =
    get_typename_wrapper(O)(_lindblad_dissipator(O.data, Id_cache), SuperOperator(), O.dimensions)

# It is already a SuperOperator
lindblad_dissipator(O::AbstractQuantumObject{SuperOperator}, Id_cache = nothing) = O

@doc raw"""
    liouvillian(H::AbstractQuantumObject, c_ops::Union{Nothing,AbstractVector,Tuple}=nothing, Id_cache=I(prod(H.dimensions)))

Construct the Liouvillian [`SuperOperator`](@ref) for a system Hamiltonian ``\hat{H}`` and a set of collapse operators ``\{\hat{C}_n\}_n``:

```math
\mathcal{L} [\cdot] = -i[\hat{H}, \cdot] + \sum_n \mathcal{D}(\hat{C}_n) [\cdot]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\cdot] = \hat{C}_n [\cdot] \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n [\cdot] - \frac{1}{2} [\cdot] \hat{C}_n^\dagger \hat{C}_n
```

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when the same function is applied multiple times with a known Hilbert space dimension.

See also [`spre`](@ref), [`spost`](@ref), and [`lindblad_dissipator`](@ref).
"""
function liouvillian(
    H::AbstractQuantumObject{OpType},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    Id_cache = I(prod(H.dimensions)),
) where {OpType<:Union{Operator,SuperOperator}}
    L = liouvillian(H, Id_cache)
    if !(c_ops isa Nothing)
        L += _sum_lindblad_dissipators(c_ops, Id_cache)
    end
    return L
end

liouvillian(H::Nothing, c_ops::Union{AbstractVector,Tuple}, Id_cache::Diagonal = I(prod(c_ops[1].dims))) =
    _sum_lindblad_dissipators(c_ops, Id_cache)

liouvillian(H::Nothing, c_ops::Nothing) = 0

liouvillian(H::AbstractQuantumObject{Operator}, Id_cache::Diagonal = I(prod(H.dimensions))) =
    get_typename_wrapper(H)(_liouvillian(H.data, Id_cache), SuperOperator(), H.dimensions)

liouvillian(H::AbstractQuantumObject{SuperOperator}, Id_cache::Diagonal) = H

function _sum_lindblad_dissipators(c_ops, Id_cache::Diagonal)
    D = 0
    # sum all the (time-independent) c_ops first
    c_ops_ti = filter(op -> isa(op, QuantumObject), c_ops)
    if !isempty(c_ops_ti)
        D += mapreduce(op -> lindblad_dissipator(op, Id_cache), +, c_ops_ti)
    end

    # sum rest of the QobjEvo together
    c_ops_td = filter(op -> isa(op, QuantumObjectEvolution), c_ops)
    if !isempty(c_ops_td)
        D += mapreduce(op -> lindblad_dissipator(op, Id_cache), +, c_ops_td)
    end

    return D
end
