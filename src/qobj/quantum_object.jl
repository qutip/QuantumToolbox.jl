#=
This file defines the QuantumObject (Qobj) structure.
It also implements the fundamental functions in Julia standard library:
    - Base: show, real, imag, Vector, Matrix
    - SparseArrays: sparse, nnz, nonzeros, rowvals, droptol!, dropzeros, dropzeros!, SparseVector, SparseMatrixCSC
    - SciMLOperators: cache_operator
=#

export QuantumObject

@doc raw"""
    struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType,N}
        data::MT
        type::ObjType
        dims::SVector{N, Int}
    end

Julia struct representing any quantum objects.

# Examples

```jldoctest
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

function Base.show(io::IO, QO::QuantumObject)
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

Base.real(x::QuantumObject) = QuantumObject(real(x.data), x.type, x.dims)
Base.imag(x::QuantumObject) = QuantumObject(imag(x.data), x.type, x.dims)

SparseArrays.sparse(A::QuantumObject{<:AbstractArray{T}}) where {T} = QuantumObject(sparse(A.data), A.type, A.dims)
SparseArrays.nnz(A::QuantumObject{<:AbstractSparseArray}) = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject{<:AbstractSparseArray}) = nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject{<:AbstractSparseArray}) = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject{<:AbstractSparseArray}, tol::Real) = (droptol!(A.data, tol); return A)
SparseArrays.dropzeros(A::QuantumObject{<:AbstractSparseArray}) = QuantumObject(dropzeros(A.data), A.type, A.dims)
SparseArrays.dropzeros!(A::QuantumObject{<:AbstractSparseArray}) = (dropzeros!(A.data); return A)

@doc raw"""
    SciMLOperators.cached_operator(L::AbstractQuantumObject, u)

Allocate caches for [`AbstractQuantumObject`](@ref) `L` for in-place evaluation with `u`-like input vectors.

Here, `u` can be in either the following types:
- `AbstractVector`
- [`Ket`](@ref)-type [`QuantumObject`](@ref) (if `L` is an [`Operator`](@ref))
- [`OperatorKet`](@ref)-type [`QuantumObject`](@ref) (if `L` is a [`SuperOperator`](@ref))
"""
SciMLOperators.cache_operator(
    L::AbstractQuantumObject{DT,OpType},
    u::AbstractVector,
) where {DT,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    get_typename_wrapper(L)(cache_operator(L.data, sparse_to_dense(similar(u))), L.type, L.dims)

function SciMLOperators.cache_operator(
    L::AbstractQuantumObject{DT1,OpType},
    u::QuantumObject{DT2,SType},
) where {
    DT1,
    DT2,
    OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    SType<:Union{KetQuantumObject,OperatorKetQuantumObject},
}
    check_dims(L, u)

    if isoper(L) && isoperket(u)
        throw(ArgumentError("The input state `u` must be a Ket if `L` is an Operator."))
    elseif issuper(L) && isket(u)
        throw(ArgumentError("The input state `u` must be an OperatorKet if `L` is a SuperOperator."))
    end
    return cache_operator(L, u.data)
end

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
