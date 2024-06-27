#=
This file defines:
    1. the QuantumObject (Qobj) structure
    2. all the type structures for QuantumObject
Also support for fundamental functions in Julia standard library:
    - Base: show, length, size, eltype, getindex, setindex!, isequal, :(==), isapprox, Vector, Matrix
    - SparseArrays: sparse, nnz, nonzeros, rowvals, droptol!, dropzeros, dropzeros!, SparseVector, SparseMatrixCSC
=#

export AbstractQuantumObject, QuantumObject
export QuantumObjectType,
    BraQuantumObject,
    KetQuantumObject,
    OperatorQuantumObject,
    OperatorBraQuantumObject,
    OperatorKetQuantumObject,
    SuperOperatorQuantumObject
export Bra, Ket, Operator, OperatorBra, OperatorKet, SuperOperator

abstract type AbstractQuantumObject end
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
    struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType}
        data::MT
        type::ObjType
        dims::Vector{Int}
    end

Julia struct representing any quantum objects.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> a isa QuantumObject
true
```
"""
struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType} <: AbstractQuantumObject
    data::MT
    type::ObjType
    dims::Vector{Int}
end

function QuantumObject(
    A::AbstractMatrix{T};
    type::ObjType = nothing,
    dims = nothing,
) where {T,ObjType<:Union{Nothing,QuantumObjectType}}
    _size = size(A)

    if type isa Nothing
        type = (_size[1] == 1 && _size[2] > 1) ? Bra : Operator # default type
    elseif type != Operator && type != SuperOperator && type != Bra && type != OperatorBra
        throw(
            ArgumentError(
                "The argument type must be Operator, SuperOperator, Bra or OperatorBra if the input array is a matrix.",
            ),
        )
    end

    if dims isa Nothing
        if type isa OperatorQuantumObject || type isa BraQuantumObject
            dims = [_size[2]]
        elseif type isa SuperOperatorQuantumObject || type isa OperatorBraQuantumObject
            dims = [isqrt(_size[2])]
        end
    else
        _check_dims(dims)
    end

    _check_QuantumObject(type, dims, _size[1], _size[2])
    return QuantumObject(A, type, dims)
end

function QuantumObject(
    A::AbstractVector{T};
    type::ObjType = nothing,
    dims = nothing,
) where {T,ObjType<:Union{Nothing,QuantumObjectType}}
    if type isa Nothing
        type = Ket # default type
    elseif type != Ket && type != OperatorKet
        throw(ArgumentError("The argument type must be Ket or OperatorKet if the input array is a vector."))
    end

    _size = (length(A), 1)

    if dims isa Nothing
        if type isa KetQuantumObject
            dims = [_size[1]]
        elseif type isa OperatorKetQuantumObject
            dims = [isqrt(_size[1])]
        end
    else
        _check_dims(dims)
    end

    _check_QuantumObject(type, dims, _size[1], _size[2])
    return QuantumObject(A, type, dims)
end

function QuantumObject(
    A::AbstractArray{T,N};
    type::ObjType = nothing,
    dims = nothing,
) where {T,N,ObjType<:Union{Nothing,QuantumObjectType}}
    throw(DomainError(size(A), "The size of the array is not compatible with vector or matrix."))
end

_check_dims(dims::Vector{Int}) =
    all(>(0), dims) || throw(DomainError(dims, "The argument dims must be a vector with positive integers."))
_check_dims(dims::Any) = throw(ArgumentError("The argument dims must be a vector with positive integers."))

function _check_QuantumObject(type::KetQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Ket"))
    (prod(dims) != m) && throw(DimensionMismatch("Ket with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::BraQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Bra"))
    (prod(dims) != n) && throw(DimensionMismatch("Bra with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with Operator"))
    (prod(dims) != m) &&
        throw(DimensionMismatch("Operator with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::SuperOperatorQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with SuperOperator"))
    (prod(dims) != sqrt(m)) &&
        throw(DimensionMismatch("SuperOperator with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorKetQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorKet"))
    (prod(dims) != sqrt(m)) &&
        throw(DimensionMismatch("OperatorKet with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function _check_QuantumObject(type::OperatorBraQuantumObject, dims::Vector{Int}, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorBra"))
    (prod(dims) != sqrt(n)) &&
        throw(DimensionMismatch("OperatorBra with dims = $(dims) does not fit the array size = $((m, n))."))
    return nothing
end

function QuantumObject(
    A::QuantumObject{<:AbstractArray{T,N}};
    type::ObjType = A.type,
    dims = A.dims,
) where {T,N,ObjType<:QuantumObjectType}
    _size = N == 1 ? (length(A), 1) : size(A)
    _check_QuantumObject(type, dims, _size[1], _size[2])
    return QuantumObject(copy(A.data), type, dims)
end

function Base.show(
    io::IO,
    QO::QuantumObject{<:AbstractArray{T},OpType},
) where {
    T,
    OpType<:Union{
        BraQuantumObject,
        KetQuantumObject,
        OperatorBraQuantumObject,
        OperatorKetQuantumObject,
        SuperOperatorQuantumObject,
    },
}
    op_data = QO.data
    println(io, "Quantum Object:   type=", QO.type, "   dims=", QO.dims, "   size=", size(op_data))
    return show(io, MIME("text/plain"), op_data)
end

function Base.show(io::IO, QO::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:OperatorQuantumObject}
    op_data = QO.data
    println(
        io,
        "Quantum Object:   type=",
        QO.type,
        "   dims=",
        QO.dims,
        "   size=",
        size(op_data),
        "   ishermitian=",
        ishermitian(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

@doc raw"""
    size(A::QuantumObject)
    size(A::QuantumObject, idx::Int)

Returns a tuple containing each dimensions of the array in the [`QuantumObject`](@ref).

Optionally, you can specify an index (`idx`) to just get the corresponding dimension of the array.
"""
Base.size(A::QuantumObject{<:AbstractArray{T}}) where {T} = size(A.data)
Base.size(A::QuantumObject{<:AbstractArray{T}}, idx::Int) where {T} = size(A.data, idx)

Base.getindex(A::QuantumObject{<:AbstractArray{T}}, inds...) where {T} = getindex(A.data, inds...)
Base.setindex!(A::QuantumObject{<:AbstractArray{T}}, val, inds...) where {T} = setindex!(A.data, val, inds...)

@doc raw"""
    eltype(A::QuantumObject)

Returns the elements type of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.eltype(A::QuantumObject) = eltype(A.data)

@doc raw"""
    length(A::QuantumObject)

Returns the length of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.length(A::QuantumObject{<:AbstractArray{T}}) where {T} = length(A.data)

Base.isequal(A::QuantumObject{<:AbstractArray{T}}, B::QuantumObject{<:AbstractArray{T}}) where {T} =
    isequal(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.isapprox(A::QuantumObject{<:AbstractArray{T}}, B::QuantumObject{<:AbstractArray{T}}) where {T} =
    isapprox(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.:(==)(A::QuantumObject{<:AbstractArray{T}}, B::QuantumObject{<:AbstractArray{T}}) where {T} =
    (A.data == B.data) && (A.type == B.type) && (A.dims == B.dims)

SparseArrays.sparse(A::QuantumObject{<:AbstractArray{T}}) where {T} = QuantumObject(sparse(A.data), A.type, A.dims)
SparseArrays.nnz(A::QuantumObject{<:AbstractSparseArray}) = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject{<:AbstractSparseArray}) = nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject{<:AbstractSparseArray}) = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject{<:AbstractSparseArray}, tol::Real) = (droptol!(A.data, tol); return A)
SparseArrays.dropzeros(A::QuantumObject{<:AbstractSparseArray}) = QuantumObject(dropzeros(A.data), A.type, A.dims)
SparseArrays.dropzeros!(A::QuantumObject{<:AbstractSparseArray}) = (dropzeros!(A.data); return A)

# data type conversions
Base.Vector(A::QuantumObject{<:AbstractVector}) = QuantumObject(Vector(A.data), A.type, A.dims)
Base.Vector{T}(A::QuantumObject{<:AbstractVector}) where {T<:Number} = QuantumObject(Vector{T}(A.data), A.type, A.dims)
Base.Matrix(A::QuantumObject{<:AbstractMatrix}) = QuantumObject(Matrix(A.data), A.type, A.dims)
Base.Matrix{T}(A::QuantumObject{<:AbstractMatrix}) where {T<:Number} = QuantumObject(Matrix{T}(A.data), A.type, A.dims)
SparseArrays.SparseVector(A::QuantumObject{<:AbstractVector}) = QuantumObject(SparseVector(A.data), A.type, A.dims)
SparseArrays.SparseVector{T}(A::QuantumObject{<:SparseVector}) where {T<:Number} =
    QuantumObject(SparseVector{T}(A.data), A.type, A.dims)
SparseArrays.SparseMatrixCSC(A::QuantumObject{<:AbstractMatrix}) =
    QuantumObject(SparseMatrixCSC(A.data), A.type, A.dims)
SparseArrays.SparseMatrixCSC{T}(A::QuantumObject{<:SparseMatrixCSC}) where {T<:Number} =
    QuantumObject(SparseMatrixCSC{T}(A.data), A.type, A.dims)
