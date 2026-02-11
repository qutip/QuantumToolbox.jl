#=  
This file defines the AbstractQuantumObject structure, all the type structures for AbstractQuantumObject, and fundamental functions in Julia standard library:  
    - Base: show, length, size, copy, eltype, getindex, setindex!, isequal, :(==), isapprox  
=#

export AbstractQuantumObject
export QuantumObjectType, SuperOperatorType, Bra, Ket, Operator, OperatorBra, OperatorKet, SuperOperator

@doc raw"""
    abstract type AbstractQuantumObject{ObjType,DimType,DataType}

Abstract type for all quantum objects like [`QuantumObject`](@ref) and [`QuantumObjectEvolution`](@ref).

# Example
```jldoctest
julia> sigmax() isa AbstractQuantumObject
true
```
"""
abstract type AbstractQuantumObject{ObjType, DimType, DataType} end

abstract type QuantumObjectType end

abstract type SuperOperatorType <: QuantumObjectType end

@doc raw"""
    Bra <: QuantumObjectType

Constructor representing a bra state ``\langle\psi|``.
"""
struct Bra <: QuantumObjectType end

Base.show(io::IO, ::Bra) = print(io, "Bra()")

@doc raw"""
    Ket <: QuantumObjectType

Constructor representing a ket state ``|\psi\rangle``.
"""
struct Ket <: QuantumObjectType end

Base.show(io::IO, ::Ket) = print(io, "Ket()")

@doc raw"""
    Operator <: QuantumObjectType

Constructor representing an operator ``\hat{O}``.
"""
struct Operator <: QuantumObjectType end

Base.show(io::IO, ::Operator) = print(io, "Operator()")

@doc raw"""
    SuperOperator <: SuperOperatorType

Constructor representing a super-operator ``\hat{\mathcal{O}}`` acting on vectorized density operator matrices.
"""
struct SuperOperator <: SuperOperatorType end

Base.show(io::IO, ::SuperOperator) = print(io, "SuperOperator()")

@doc raw"""
    OperatorBra <: QuantumObjectType

Constructor representing a bra state in the [`SuperOperator`](@ref) formalism ``\langle\!\langle\rho|``.
"""
struct OperatorBra <: QuantumObjectType end

Base.show(io::IO, ::OperatorBra) = print(io, "OperatorBra()")

@doc raw"""
    OperatorKet <: QuantumObjectType

Constructor representing a ket state in the [`SuperOperator`](@ref) formalism ``|\rho\rangle\!\rangle``.
"""
struct OperatorKet <: QuantumObjectType end

Base.show(io::IO, ::OperatorKet) = print(io, "OperatorKet()")

@doc raw"""
    size(A::AbstractQuantumObject)
    size(A::AbstractQuantumObject, idx::Int)
    shape(A::AbstractQuantumObject)
    shape(A::AbstractQuantumObject, idx::Int)

Returns a tuple containing each dimensions of the array in the [`AbstractQuantumObject`](@ref).

Optionally, you can specify an index (`idx`) to just get the corresponding dimension of the array.

!!! note
    `shape` is a synonym of `size`.
"""
Base.size(A::AbstractQuantumObject) = size(A.data)
Base.size(A::AbstractQuantumObject, idx::Int) = size(A.data, idx)

Base.copy(A::AbstractQuantumObject) = get_typename_wrapper(A)(copy(A.data), A.type, A.dimensions)

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
    isequal(A.type, B.type) && isequal(A.dimensions, B.dimensions) && isequal(A.data, B.data)
Base.isapprox(A::AbstractQuantumObject, B::AbstractQuantumObject; kwargs...) =
    isequal(A.type, B.type) && isequal(A.dimensions, B.dimensions) && isapprox(A.data, B.data; kwargs...)
Base.:(==)(A::AbstractQuantumObject, B::AbstractQuantumObject) =
    (A.type == B.type) && (A.dimensions == B.dimensions) && (A.data == B.data)

function check_dimensions(dimensions_list::NTuple{N, Dimensions}) where {N}
    allequal(dimensions_list) ||
        throw(DimensionMismatch("The quantum objects should have the same Hilbert `dimensions`."))
    return nothing
end

function check_dimensions(Qobj_tuple::NTuple{N, AbstractQuantumObject}) where {N}
    dimensions_list = map(A -> A.dimensions, Qobj_tuple)
    allequal(dimensions_list) ||
        throw(DimensionMismatch("The quantum objects should have the same Hilbert `dimensions`."))
    return nothing
end
check_dimensions(A::AbstractQuantumObject...) = check_dimensions(A)

function _check_QuantumObject(::Ket, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[2] != 1) && throw(DimensionMismatch(("The size $(array_size) of the array is not compatible with Ket")))
    _check_dims_and_array_size("Ket", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::Bra, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[1] != 1) && throw(DimensionMismatch(("The size $(array_size) of the array is not compatible with Bra")))
    _check_dims_and_array_size("Bra", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::Operator, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    _check_dims_and_array_size("Operator", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::SuperOperator, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    _check_dims_and_array_size("SuperOperator", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::OperatorKet, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[2] != 1) && throw(DimensionMismatch(("The size $(array_size) of the array is not compatible with OperatorKet")))
    _check_dims_and_array_size("OperatorKet", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::OperatorBra, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[1] != 1) && throw(DimensionMismatch(("The size $(array_size) of the array is not compatible with OperatorBra")))
    _check_dims_and_array_size("OperatorBra", dimensions, array_size)
    return nothing
end

# this helps to check if the Dimensions matches the array size
function _check_dims_and_array_size(type_string::String, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (get_size(dimensions) == array_size) || throw(
        DimensionMismatch(
            "$(type_string) with dims = $(_get_dims_string(dimensions)) does not fit the array size = $(array_size).",
        ),
    )
    return nothing
end

_check_type(::T) where {T <: Union{Nothing, <:QuantumObjectType}} = T
_check_type(::Type{T}) where {T} =
    throw(ArgumentError("The argument `$T` is not valid. You may probably want to use `$T()` instead."))
_check_type(t) = throw(ArgumentError("The argument $t is not valid. It should be a subtype of `QuantumObjectType`."))

function Base.getproperty(A::AbstractQuantumObject, key::Symbol)
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(A, :dimensions))
    else
        return getfield(A, key)
    end
end

# _gen_data_size help us handle the data size
# especially for AbstractVector, it returns a 2-element tuple, which is necessary in _check_dims_and_array_size
_gen_data_size(data::AbstractVector) = (size(data, 1), 1)
_gen_data_size(data::Union{AbstractArray, AbstractSciMLOperator}) = size(data)

# generate dimensions based on different QuantumObjectType and different input formats of dims:
## dims::Integer
## (backward compatibility, but avoid support for OperatorKet/OperatorBra/SuperOperator since it causes ambiguity)
_gen_dimensions(::Ket, dims::Integer) = Dimensions(Space(dims), Space(1))
_gen_dimensions(::Bra, dims::Integer) = Dimensions(Space(1), Space(dims))
_gen_dimensions(::Operator, dims::Integer) = Dimensions(Space(dims)) # is endomorphism

## dims::VectorOrTuple{Int} : vector or tuple of integers
## (backward compatibility, but avoid support for OperatorKet/OperatorBra/SuperOperator since it causes ambiguity)
_gen_dimensions(::Ket, dims::VectorOrTuple{T}) where {T <: Integer} = Dimensions(_list_to_tensor_space(dims, "dims"), Space(1))
_gen_dimensions(::Bra, dims::VectorOrTuple{T}) where {T <: Integer} = Dimensions(Space(1), _list_to_tensor_space(dims, "dims"))
_gen_dimensions(::Operator, dims::VectorOrTuple{T}) where {T <: Integer} = Dimensions(_list_to_tensor_space(dims, "dims"))

## dims::DimsListType{T1, T2} : general cases
_gen_dimensions(::QuantumObjectType, dims::DimsListType{T1, T2}) where {T1, T2} = Dimensions(dims)

# other cases
_gen_dimensions(::QuantumObjectType, dims::Dimensions) = dims
_gen_dimensions(type, dims) = throw(ArgumentError("The argument `dims` with value $dims is not valid for object type $type."))

# function _gen_dimensions(::ObjType, dims::Union{T, VectorOrTuple{T}}) where {T <: Integer}

#     if ObjType <: Union{Ket, OperatorKet}
#         return Dimensions(raw_dimensions.to, Space(1))
#     elseif ObjType <: Union{Bra, OperatorBra}
#         return Dimensions(Space(1), raw_dimensions.from)
#     else
#         return raw_dimensions
#     end
# end
# _gen_dimensions(type::QuantumObjectType, dims::Union{AbstractVector{T}, NTuple{N, T}}) where {T <: Union{AbstractVector, NTuple}, N} =
#     Dimensions(dims)
# _gen_dimensions(type, dims) = throw(ArgumentError("The argument `dims` with value $dims is not valid for object type $type."))

# functions for getting Float or Complex element type
_float_type(A::AbstractQuantumObject) = _float_type(eltype(A))
_complex_float_type(A::AbstractQuantumObject) = _complex_float_type(eltype(A))
