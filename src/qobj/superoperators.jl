#=
Functions for generating (common) quantum super-operators.
=#

export spre, spost, sprepost, liouvillian, lindblad_dissipator

# intrinsic functions for super-operators
## keep these because they take AbstractMatrix as input and ensure the output is sparse matrix
_spre(A::AbstractMatrix) = kron(Eye(size(A, 1)), sparse(A))
_spre(A::AbstractSparseMatrix) = kron(Eye(size(A, 1)), A)
_spost(B::AbstractMatrix) = kron(transpose(sparse(B)), Eye(size(B, 1)))
_spost(B::AbstractSparseMatrix) = kron(transpose(B), Eye(size(B, 1)))
_sprepost(A::AbstractMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), sparse(A))
_sprepost(A::AbstractMatrix, B::AbstractSparseMatrix) = kron(transpose(B), sparse(A))
_sprepost(A::AbstractSparseMatrix, B::AbstractMatrix) = kron(transpose(sparse(B)), A)
_sprepost(A::AbstractSparseMatrix, B::AbstractSparseMatrix) = kron(transpose(B), A)
_sprepost(A, B) = _spre(A) * _spost(B) # for any other input types

## if input is AbstractSciMLOperator
## some of them are optimized to speed things up
## the rest of the SciMLOperators will just use lazy tensor (and prompt a warning)
_spre(A::MatrixOperator) = MatrixOperator(_spre(A.A))
_spre(A::ScaledOperator) = ScaledOperator(A.λ, _spre(A.L))
_spre(A::AddedOperator) = AddedOperator(map(op -> _spre(op), A.ops))
_spre(A::ComposedOperator) = ComposedOperator(map(_spre, A.ops), nothing)
function _spre(A::AbstractSciMLOperator)
    Id = Eye(size(A, 1))
    _lazy_tensor_warning(Id, A)
    return kron(Id, A)
end

_spost(B::MatrixOperator) = MatrixOperator(_spost(B.A))
_spost(B::ScaledOperator) = ScaledOperator(B.λ, _spost(B.L))
_spost(B::AddedOperator) = AddedOperator(map(op -> _spost(op), B.ops))
_spost(B::ComposedOperator) = ComposedOperator(map(_spost, reverse(B.ops)), nothing)
function _spost(B::AbstractSciMLOperator)
    B_T = transpose(B)
    Id = Eye(size(B, 1))
    _lazy_tensor_warning(B_T, Id)
    return kron(B_T, Id)
end

## intrinsic liouvillian
function _liouvillian(
        H::MT,
        ::Val{assume_hermitian},
    ) where {MT <: Union{AbstractMatrix, AbstractSciMLOperator}, assume_hermitian}
    H_spre = _spre(H)
    H_spost = assume_hermitian ? _spost(H) : _spost(H')
    return -1im * (H_spre - H_spost)
end
function _liouvillian(H::MatrixOperator, assume_hermitian::Val)
    isconstant(H) ||
        throw(ArgumentError("The given Hamiltonian for constructing Liouvillian must be constant MatrixOperator."))
    return MatrixOperator(_liouvillian(H.A, assume_hermitian))
end
_liouvillian(H::ScaledOperator, assume_hermitian::Val{true}) = ScaledOperator(H.λ, _liouvillian(H.L, assume_hermitian))
_liouvillian(H::AddedOperator, assume_hermitian::Val) =
    AddedOperator(map(op -> _liouvillian(op, assume_hermitian), H.ops))

# intrinsic lindblad_dissipator
function _lindblad_dissipator(O::MT) where {MT <: Union{AbstractMatrix, AbstractSciMLOperator}}
    Od_O = O' * O
    return _sprepost(O, O') - (_spre(Od_O) + _spost(Od_O)) / 2
end
function _lindblad_dissipator(O::MatrixOperator)
    _O = O.A
    Od_O = _O' * _O
    return MatrixOperator(_sprepost(_O, _O') - (_spre(Od_O) + _spost(Od_O)) / 2)
end
function _lindblad_dissipator(O::ScaledOperator)
    λc_λ = conj(O.λ) * O.λ
    return ScaledOperator(λc_λ, _lindblad_dissipator(O.L))
end

@doc raw"""
    spre(A::AbstractQuantumObject)

Returns the [`SuperOperator`](@ref) form of `A` acting on the left of the density matrix operator: ``\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\!\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{\mathbb{1}} \otimes \hat{A}``, namely 

```math
\mathcal{O} \left(\hat{A}\right) \left[ \hat{\rho} \right] = \hat{\mathbb{1}} \otimes \hat{A} ~ |\hat{\rho}\rangle\!\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

See also [`spost`](@ref) and [`sprepost`](@ref).
"""
spre(A::AbstractQuantumObject{Operator}) = get_typename_wrapper(A)(_spre(A.data), SuperOperator(), A.dimensions)

@doc raw"""
    spost(B::AbstractQuantumObject)

Returns the [`SuperOperator`](@ref) form of `B` acting on the right of the density matrix operator: ``\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\!\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{\mathbb{1}}``, namely

```math
\mathcal{O} \left(\hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{\mathbb{1}} ~ |\hat{\rho}\rangle\!\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

See also [`spre`](@ref) and [`sprepost`](@ref).
"""
spost(B::AbstractQuantumObject{Operator}) = get_typename_wrapper(B)(_spost(B.data), SuperOperator(), B.dimensions)

@doc raw"""
    sprepost(A::AbstractQuantumObject, B::AbstractQuantumObject)

Returns the [`SuperOperator`](@ref) form of `A` and `B` acting on the left and right of the density matrix operator, respectively: ``\mathcal{O} \left( \hat{A}, \hat{B} \right) \left[ \hat{\rho} \right] = \hat{A} \hat{\rho} \hat{B}``.

Since the density matrix is vectorized in [`OperatorKet`](@ref) form: ``|\hat{\rho}\rangle\!\rangle``, this [`SuperOperator`](@ref) is always a matrix ``\hat{B}^T \otimes \hat{A}``, namely

```math
\mathcal{O} \left(\hat{A}, \hat{B}\right) \left[ \hat{\rho} \right] = \hat{B}^T \otimes \hat{A} ~ |\hat{\rho}\rangle\!\rangle = \textrm{spre}(\hat{A}) * \textrm{spost}(\hat{B}) ~ |\hat{\rho}\rangle\!\rangle
```
(see the section in documentation: [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators) for more details)

See also [`spre`](@ref) and [`spost`](@ref).
"""
function sprepost(A::AbstractQuantumObject{Operator}, B::AbstractQuantumObject{Operator})
    check_dimensions(A, B)
    return promote_op_type(A, B)(_sprepost(A.data, B.data), SuperOperator(), A.dimensions)
end

@doc raw"""
    lindblad_dissipator(O::AbstractQuantumObject)

Returns the Lindblad [`SuperOperator`](@ref) defined as

```math
\mathcal{D} \left( \hat{O} \right) \left[ \hat{\rho} \right] = \frac{1}{2} \left( 2 \hat{O} \hat{\rho} \hat{O}^\dagger - 
\hat{O}^\dagger \hat{O} \hat{\rho} - \hat{\rho} \hat{O}^\dagger \hat{O} \right)
```

See also [`spre`](@ref), [`spost`](@ref), and [`sprepost`](@ref).
"""
lindblad_dissipator(O::AbstractQuantumObject{Operator}) =
    get_typename_wrapper(O)(_lindblad_dissipator(O.data), SuperOperator(), O.dimensions)

# It is already a SuperOperator
lindblad_dissipator(O::AbstractQuantumObject{SuperOperator}) = O

@doc raw"""
    liouvillian(
        H::AbstractQuantumObject,
        c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        assume_hermitian::Union{Bool,Val} = Val(true),
    )

Construct the Liouvillian [`SuperOperator`](@ref) for a system Hamiltonian ``\hat{H}`` and a set of collapse operators ``\{\hat{C}_n\}_n``.

By default, when the Hamiltonian `H` is assumed to be Hermitian [`assume_hermitian = Val(true)` or `true`], the Liouvillian [`SuperOperator`](@ref) is defined as :

```math
\mathcal{L} [\cdot] = -i\left(\hat{H}[\cdot] - [\cdot]\hat{H}\right) + \sum_n \mathcal{D}(\hat{C}_n) [\cdot],
```

otherwise [`assume_hermitian = Val(false)` or `false`],

```math
\mathcal{L} [\cdot] = -i\left(\hat{H}[\cdot] - [\cdot]\hat{H}^\dagger\right) + \sum_n \mathcal{D}(\hat{C}_n) [\cdot],
```

where 

```math
\mathcal{D}(\hat{C}_n) [\cdot] = \hat{C}_n [\cdot] \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n [\cdot] - \frac{1}{2} [\cdot] \hat{C}_n^\dagger \hat{C}_n
```

See also [`spre`](@ref), [`spost`](@ref), and [`lindblad_dissipator`](@ref).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `assume_hermitian = Val(true)` instead of `assume_hermitian = true`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function liouvillian(
        H::AbstractQuantumObject{OpType},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        assume_hermitian::Union{Bool, Val} = Val(true),
    ) where {OpType <: Union{Operator, SuperOperator}}
    L = liouvillian(H; assume_hermitian = assume_hermitian)
    if !isnothing(c_ops)
        return L + _sum_lindblad_dissipators(c_ops)
    end
    return L
end

liouvillian(H::Nothing, c_ops::Union{AbstractVector, Tuple}; kwargs...) = _sum_lindblad_dissipators(c_ops)

liouvillian(H::Nothing, c_ops::Nothing; kwargs...) = 0

liouvillian(H::AbstractQuantumObject{Operator}; assume_hermitian::Union{Bool, Val} = Val(true)) =
    get_typename_wrapper(H)(_liouvillian(H.data, makeVal(assume_hermitian)), SuperOperator(), H.dimensions)

liouvillian(H::AbstractQuantumObject{SuperOperator}; kwargs...) = H
_sum_lindblad_dissipators(c_ops::Nothing) = 0

_sum_lindblad_dissipators(c_ops::AbstractVector) = sum(op -> lindblad_dissipator(op), c_ops; init = 0)

# Help the compiler to unroll the sum at compile time
@generated function _sum_lindblad_dissipators(c_ops::Tuple)
    N = length(c_ops.parameters)
    if N == 0
        return :(0)
    end
    ex = :(lindblad_dissipator(c_ops[1]))
    for i in 2:N
        ex = :($ex + lindblad_dissipator(c_ops[$i]))
    end
    return ex
end
