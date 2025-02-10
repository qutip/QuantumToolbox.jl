#=
This file defines the QuantumObject (Qobj) structure.
It also implements the fundamental functions in Julia standard library:
    - Base: show, real, imag, Vector, Matrix
    - SparseArrays: sparse, nnz, nonzeros, rowvals, droptol!, dropzeros, dropzeros!, SparseVector, SparseMatrixCSC
    - SciMLOperators: cache_operator
=#

export QuantumObject

@doc raw"""
    struct QuantumObject{ObjType<:QuantumObjectType,DimType<:AbstractDimensions,DataType<:AbstractArray} <: AbstractQuantumObject{ObjType,DimType,DataType}
        data::DataType
        type::ObjType
        dimensions::DimType
    end

Julia structure representing any time-independent quantum objects. For time-dependent cases, see [`QuantumObjectEvolution`](@ref).

!!! note "`dims` property"
    For a given `H::QuantumObject`, `H.dims` or `getproperty(H, :dims)` returns its `dimensions` in the type of integer-vector.

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

julia> a.dims
1-element SVector{1, Int64} with indices SOneTo(1):
 20

julia> a.dimensions
Dimensions{1, Tuple{Space}}((Space(20),))
```
"""
struct QuantumObject{ObjType<:QuantumObjectType,DimType<:AbstractDimensions,DataType<:AbstractArray} <:
       AbstractQuantumObject{ObjType,DimType,DataType}
    data::DataType
    type::ObjType
    dimensions::DimType

    function QuantumObject(data::DT, type::ObjType, dims) where {DT<:AbstractArray,ObjType<:QuantumObjectType}
        dimensions = _gen_dimensions(dims)

        _size = _get_size(data)
        _check_QuantumObject(type, dimensions, _size[1], _size[2])

        return new{ObjType,typeof(dimensions),DT}(data, type, dimensions)
    end
end

@doc raw"""
    Qobj(A::AbstractArray; type = nothing, dims = nothing)
    QuantumObject(A::AbstractArray; type = nothing, dims = nothing)

Generate [`QuantumObject`](@ref) with a given `A::AbstractArray` and specified `type::QuantumObjectType` and `dims`.

!!! note
    `Qobj` is a synonym of `QuantumObject`.
"""
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
        if type isa BraQuantumObject
            dims = Dimensions(_size[2])
        elseif type isa OperatorQuantumObject
            dims =
                (_size[1] == _size[2]) ? Dimensions(_size[1]) :
                GeneralDimensions(SVector{2}(SVector{1}(_size[1]), SVector{1}(_size[2])))
        elseif type isa SuperOperatorQuantumObject || type isa OperatorBraQuantumObject
            dims = Dimensions(isqrt(_size[2]))
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
            dims = Dimensions(_size[1])
        elseif type isa OperatorKetQuantumObject
            dims = Dimensions(isqrt(_size[1]))
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

function QuantumObject(A::QuantumObject; type::ObjType = A.type, dims = A.dimensions) where {ObjType<:QuantumObjectType}
    _size = _get_size(A.data)
    dimensions = _gen_dimensions(dims)
    _check_QuantumObject(type, dimensions, _size[1], _size[2])
    return QuantumObject(copy(A.data), type, dimensions)
end

function Base.show(
    io::IO,
    QO::QuantumObject{OpType},
) where {
    OpType<:Union{
        BraQuantumObject,
        KetQuantumObject,
        OperatorBraQuantumObject,
        OperatorKetQuantumObject,
        SuperOperatorQuantumObject,
    },
}
    op_data = QO.data
    println(
        io,
        "\nQuantum Object:   type=",
        QO.type,
        "   dims=",
        _get_dims_string(QO.dimensions),
        "   size=",
        size(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

function Base.show(io::IO, QO::QuantumObject)
    op_data = QO.data
    println(
        io,
        "\nQuantum Object:   type=",
        QO.type,
        "   dims=",
        _get_dims_string(QO.dimensions),
        "   size=",
        size(op_data),
        "   ishermitian=",
        ishermitian(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

Base.real(x::QuantumObject) = QuantumObject(real(x.data), x.type, x.dimensions)
Base.imag(x::QuantumObject) = QuantumObject(imag(x.data), x.type, x.dimensions)

SparseArrays.sparse(A::QuantumObject) = QuantumObject(sparse(A.data), A.type, A.dimensions)
SparseArrays.nnz(A::QuantumObject) = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject) = nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject) = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject, tol::Real) = (droptol!(A.data, tol); return A)
SparseArrays.dropzeros(A::QuantumObject) = QuantumObject(dropzeros(A.data), A.type, A.dimensions)
SparseArrays.dropzeros!(A::QuantumObject) = (dropzeros!(A.data); return A)

@doc raw"""
    SciMLOperators.cached_operator(L::AbstractQuantumObject, u)

Allocate caches for [`AbstractQuantumObject`](@ref) `L` for in-place evaluation with `u`-like input vectors.

Here, `u` can be in either the following types:
- `AbstractVector`
- [`Ket`](@ref)-type [`QuantumObject`](@ref) (if `L` is an [`Operator`](@ref))
- [`OperatorKet`](@ref)-type [`QuantumObject`](@ref) (if `L` is a [`SuperOperator`](@ref))
"""
SciMLOperators.cache_operator(
    L::AbstractQuantumObject{OpType},
    u::AbstractVector,
) where {OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    get_typename_wrapper(L)(cache_operator(L.data, to_dense(similar(u))), L.type, L.dimensions)

function SciMLOperators.cache_operator(
    L::AbstractQuantumObject{OpType},
    u::QuantumObject{SType},
) where {
    OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    SType<:Union{KetQuantumObject,OperatorKetQuantumObject},
}
    check_dimensions(L, u)

    if isoper(L) && isoperket(u)
        throw(ArgumentError("The input state `u` must be a Ket if `L` is an Operator."))
    elseif issuper(L) && isket(u)
        throw(ArgumentError("The input state `u` must be an OperatorKet if `L` is a SuperOperator."))
    end
    return cache_operator(L, u.data)
end

# data type conversions
Base.Vector(A::QuantumObject) = QuantumObject(Vector(A.data), A.type, A.dimensions)
Base.Vector{T}(A::QuantumObject) where {T<:Number} = QuantumObject(Vector{T}(A.data), A.type, A.dimensions)
Base.Matrix(A::QuantumObject) = QuantumObject(Matrix(A.data), A.type, A.dimensions)
Base.Matrix{T}(A::QuantumObject) where {T<:Number} = QuantumObject(Matrix{T}(A.data), A.type, A.dimensions)
SparseArrays.SparseVector(A::QuantumObject) = QuantumObject(SparseVector(A.data), A.type, A.dimensions)
SparseArrays.SparseVector{T}(A::QuantumObject) where {T<:Number} =
    QuantumObject(SparseVector{T}(A.data), A.type, A.dimensions)
SparseArrays.SparseMatrixCSC(A::QuantumObject) = QuantumObject(SparseMatrixCSC(A.data), A.type, A.dimensions)
SparseArrays.SparseMatrixCSC{T}(A::QuantumObject) where {T<:Number} =
    QuantumObject(SparseMatrixCSC{T}(A.data), A.type, A.dimensions)
