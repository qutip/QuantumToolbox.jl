#=  
This file defines the AbstractQuantumObject structure, all the type structures for AbstractQuantumObject, and fundamental functions in Julia standard library:  
    - Base: show, length, size, eltype, getindex, setindex!, isequal, :(==), isapprox  
=#

export AbstractQuantumObject
export QuantumObjectType,
    BraQuantumObject,
    KetQuantumObject,
    OperatorQuantumObject,
    OperatorBraQuantumObject,
    OperatorKetQuantumObject,
    SuperOperatorQuantumObject
export Bra, Ket, Operator, OperatorBra, OperatorKet, SuperOperator

@doc raw"""
    abstract type AbstractQuantumObject{DataType,ObjType,N}

Abstract type for all quantum objects like [`QuantumObject`](@ref) and [`QuantumObjectEvolution`](@ref).

# Example
```jldoctest
julia> sigmax() isa AbstractQuantumObject
true
```
"""
abstract type AbstractQuantumObject{DataType,ObjType,N} end

abstract type QuantumObjectType end

@doc raw"""
    BraQuantumObject <: QuantumObjectType

Constructor representing a bra state ``\langle\psi|``.
"""
struct BraQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::BraQuantumObject) = print(io, "Bra")

@doc raw"""
    const Bra = BraQuantumObject()

A constant representing the type of [`BraQuantumObject`](@ref): a bra state ``\langle\psi|``
"""
const Bra = BraQuantumObject()

@doc raw"""
    KetQuantumObject <: QuantumObjectType

Constructor representing a ket state ``|\psi\rangle``.
"""
struct KetQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::KetQuantumObject) = print(io, "Ket")

@doc raw"""
    const Ket = KetQuantumObject()

A constant representing the type of [`KetQuantumObject`](@ref): a ket state ``|\psi\rangle``
"""
const Ket = KetQuantumObject()

@doc raw"""
    OperatorQuantumObject <: QuantumObjectType

Constructor representing an operator ``\hat{O}``.
"""
struct OperatorQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::OperatorQuantumObject) = print(io, "Operator")

@doc raw"""
    const Operator = OperatorQuantumObject()

A constant representing the type of [`OperatorQuantumObject`](@ref): an operator ``\hat{O}``
"""
const Operator = OperatorQuantumObject()

@doc raw"""
    SuperOperatorQuantumObject <: QuantumObjectType

Constructor representing a super-operator ``\hat{\mathcal{O}}`` acting on vectorized density operator matrices.
"""
struct SuperOperatorQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::SuperOperatorQuantumObject) = print(io, "SuperOperator")

@doc raw"""
    const SuperOperator = SuperOperatorQuantumObject()

A constant representing the type of [`SuperOperatorQuantumObject`](@ref): a super-operator ``\hat{\mathcal{O}}`` acting on vectorized density operator matrices
"""
const SuperOperator = SuperOperatorQuantumObject()

@doc raw"""
    OperatorBraQuantumObject <: QuantumObjectType

Constructor representing a bra state in the [`SuperOperator`](@ref) formalism ``\langle\langle\rho|``.
"""
struct OperatorBraQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::OperatorBraQuantumObject) = print(io, "OperatorBra")

@doc raw"""
    const OperatorBra = OperatorBraQuantumObject()

A constant representing the type of [`OperatorBraQuantumObject`](@ref): a bra state in the [`SuperOperator`](@ref) formalism ``\langle\langle\rho|``.
"""
const OperatorBra = OperatorBraQuantumObject()

@doc raw"""
    OperatorKetQuantumObject <: QuantumObjectType

Constructor representing a ket state in the [`SuperOperator`](@ref) formalism ``|\rho\rangle\rangle``.
"""
struct OperatorKetQuantumObject <: QuantumObjectType end
Base.show(io::IO, ::OperatorKetQuantumObject) = print(io, "OperatorKet")

@doc raw"""
    const OperatorKet = OperatorKetQuantumObject()

A constant representing the type of [`OperatorKetQuantumObject`](@ref): a ket state in the [`SuperOperator`](@ref) formalism ``|\rho\rangle\rangle``
"""
const OperatorKet = OperatorKetQuantumObject()

@doc raw"""
    size(A::AbstractQuantumObject)
    size(A::AbstractQuantumObject, idx::Int)

Returns a tuple containing each dimensions of the array in the [`AbstractQuantumObject`](@ref).

Optionally, you can specify an index (`idx`) to just get the corresponding dimension of the array.
"""
Base.size(A::AbstractQuantumObject) = size(A.data)
Base.size(A::AbstractQuantumObject, idx::Int) = size(A.data, idx)

Base.getindex(A::AbstractQuantumObject, inds...) = getindex(A.data, inds...)
Base.setindex!(A::AbstractQuantumObject, val, inds...) = setindex!(A.data, val, inds...)

@doc raw"""
    eltype(A::AbstractQuantumObject)

Returns the elements type of the matrix or vector corresponding to the [`AbstractQuantumObject`](@ref) `A`.
"""
Base.eltype(A::AbstractQuantumObject) = eltype(A.data)

@doc raw"""
    length(A::AbstractQuantumObject)

Returns the length of the matrix or vector corresponding to the [`AbstractQuantumObject`](@ref) `A`.
"""
Base.length(A::AbstractQuantumObject) = length(A.data)

Base.isequal(A::AbstractQuantumObject, B::AbstractQuantumObject) =
    isequal(A.type, B.type) && isequal(A.dims, B.dims) && isequal(A.data, B.data)
Base.isapprox(A::AbstractQuantumObject, B::AbstractQuantumObject; kwargs...) =
    isequal(A.type, B.type) && isequal(A.dims, B.dims) && isapprox(A.data, B.data; kwargs...)
Base.:(==)(A::AbstractQuantumObject, B::AbstractQuantumObject) =
    (A.type == B.type) && (A.dims == B.dims) && (A.data == B.data)

function check_dims(A::AbstractQuantumObject, B::AbstractQuantumObject)
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects don't have the same Hilbert dimension."))
    return nothing
end

function _check_dims(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N}
    _non_static_array_warning("dims", dims)
    return (all(>(0), dims) && length(dims) > 0) ||
           throw(DomainError(dims, "The argument dims must be of non-zero length and contain only positive integers."))
end
_check_dims(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

function _check_QuantumObject(type::KetQuantumObject, dims, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Ket"))
    (prod(dims) != m) && throw(DimensionMismatch("Ket with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::BraQuantumObject, dims, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Bra"))
    (prod(dims) != n) && throw(DimensionMismatch("Bra with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorQuantumObject, dims, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with Operator"))
    (prod(dims) != m) &&
        throw(DimensionMismatch("Operator with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::SuperOperatorQuantumObject, dims, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with SuperOperator"))
    (prod(dims) != sqrt(m)) &&
        throw(DimensionMismatch("SuperOperator with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorKetQuantumObject, dims, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorKet"))
    (prod(dims) != sqrt(m)) &&
        throw(DimensionMismatch("OperatorKet with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorBraQuantumObject, dims, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorBra"))
    (prod(dims) != sqrt(n)) &&
        throw(DimensionMismatch("OperatorBra with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

# functions for getting Float or Complex element type
_FType(A::AbstractQuantumObject) = _FType(eltype(A))
_CType(A::AbstractQuantumObject) = _CType(eltype(A))
