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

@doc raw"""
    AbstractQuantumObject{DataType,ObjType,N}

Abstract type for all quantum objects like [`QuantumObject`](@ref) and [`QuantumObjectEvolution`](@ref).
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
    struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType,N}
        data::MT
        type::ObjType
        dims::SVector{N, Int}
    end

Julia struct representing any quantum objects.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⎡⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⎦

julia> a isa QuantumObject
true
```
"""
struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType,N} <: AbstractQuantumObject{MT,ObjType,N}
    data::MT
    type::ObjType
    dims::SVector{N,Int}

    function QuantumObject(data::MT, type::ObjType, dims) where {MT<:AbstractArray,ObjType<:QuantumObjectType}
        _check_dims(dims)

        _size = _get_size(data)
        _check_QuantumObject(type, dims, _size[1], _size[2])

        N = length(dims)

        return new{MT,ObjType,N}(data, type, SVector{N,Int}(dims))
    end
end

function QuantumObject(A::AbstractArray, type::ObjType, dims::Integer) where {ObjType<:QuantumObjectType}
    return QuantumObject(A, type, SVector{1,Int}(dims))
end

function QuantumObject(
    A::AbstractMatrix{T};
    type::ObjType = nothing,
    dims = nothing,
) where {T,ObjType<:Union{Nothing,QuantumObjectType}}
    _size = _get_size(A)

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
            dims = SVector{1,Int}(_size[2])
        elseif type isa SuperOperatorQuantumObject || type isa OperatorBraQuantumObject
            dims = SVector{1,Int}(isqrt(_size[2]))
        end
    end

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

    if dims isa Nothing
        _size = _get_size(A)
        if type isa KetQuantumObject
            dims = SVector{1,Int}(_size[1])
        elseif type isa OperatorKetQuantumObject
            dims = SVector{1,Int}(isqrt(_size[1]))
        end
    end

    return QuantumObject(A, type, dims)
end

function QuantumObject(
    A::AbstractArray{T,N};
    type::ObjType = nothing,
    dims = nothing,
) where {T,N,ObjType<:Union{Nothing,QuantumObjectType}}
    throw(DomainError(size(A), "The size of the array is not compatible with vector or matrix."))
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

function Base.show(io::IO, QO::AbstractQuantumObject)
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

Base.real(x::QuantumObject) = QuantumObject(real(x.data), x.type, x.dims)
Base.imag(x::QuantumObject) = QuantumObject(imag(x.data), x.type, x.dims)

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

get_typename_wrapper(A::AbstractQuantumObject) = Base.typename(typeof(A)).wrapper

# functions for getting Float or Complex element type
_FType(::QuantumObject{<:AbstractArray{T}}) where {T<:Number} = _FType(T)
_CType(::QuantumObject{<:AbstractArray{T}}) where {T<:Number} = _CType(T)
