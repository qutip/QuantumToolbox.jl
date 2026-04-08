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

function check_dimensions(Qobj_tuple::NTuple{N, AbstractQuantumObject}) where {N}
    dimensions_list = map(A -> A.dimensions, Qobj_tuple)
    allequal(dimensions_list) ||
        throw(DimensionMismatch("The quantum objects should have the same Hilbert `dimensions`."))
    return nothing
end
check_dimensions(A::AbstractQuantumObject...) = check_dimensions(A)

function _check_QuantumObject(::Ket, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (length(array_size) == 1) || _type_and_array_mismatch_error("Ket", array_size)

    # since the length of array_size is confirmed to be 1 above, here we need to make sure:
    #   - get_size(dimensions)[1] == array_size[1]
    #   - get_size(dimensions)[2] == 1
    ((get_size(dimensions)[1] == array_size[1]) && (get_size(dimensions)[2] == 1)) || _dims_and_array_mismatch_error("Ket", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::Bra, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[1] == 1) || _type_and_array_mismatch_error("Bra", array_size)
    (get_size(dimensions) == array_size) || _dims_and_array_mismatch_error("Bra", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::Operator, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (get_size(dimensions) == array_size) || _dims_and_array_mismatch_error("Operator", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::SuperOperator, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (get_size(dimensions) == array_size) || _dims_and_array_mismatch_error("SuperOperator", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::OperatorKet, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (length(array_size) == 1) || _type_and_array_mismatch_error("OperatorKet", array_size)

    # since the length of array_size is confirmed to be 1 above, here we need to make sure:
    #   - get_size(dimensions)[1] == array_size[1]
    #   - get_size(dimensions)[2] == 1
    ((get_size(dimensions)[1] == array_size[1]) && (get_size(dimensions)[2] == 1)) || _dims_and_array_mismatch_error("OperatorKet", dimensions, array_size)
    return nothing
end

function _check_QuantumObject(::OperatorBra, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N}
    (array_size[1] == 1) || _type_and_array_mismatch_error("OperatorBra", array_size)
    (get_size(dimensions) == array_size) || _dims_and_array_mismatch_error("OperatorBra", dimensions, array_size)
    return nothing
end

# these functions help to check if the type, dimensions, and array size all matches
_type_and_array_mismatch_error(type_string::String, array_size::NTuple{N, Int}) where {N} =
    throw(DimensionMismatch(("The size $(array_size) of the array is not compatible with $(type_string)")))
_dims_and_array_mismatch_error(type_string::String, dimensions::Dimensions, array_size::NTuple{N, Int}) where {N} =
    throw(
    DimensionMismatch(
        "$(type_string) with dims = $(_get_dims_string(dimensions)) does not fit the array size = $(array_size).",
    ),
)

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

# generate dimensions based on different QuantumObjectType and different input formats of dims:
## dims::Integer
## (backward compatibility, but avoid support for OperatorKet/OperatorBra/SuperOperator since it causes ambiguity)
_gen_dimensions(type::ObjType, dims::Integer) where {ObjType <: Union{Ket, Bra, Operator}} = _gen_dimensions(type, Space(dims))

## dims::AbstractVecOrTuple{Int} : vector or tuple of integers
## (backward compatibility, but avoid support for OperatorKet/OperatorBra/SuperOperator since it causes ambiguity)
_gen_dimensions(type::ObjType, dims::AbstractVecOrTuple{T}) where {ObjType <: Union{Ket, Bra, Operator}, T <: Integer} = _gen_dimensions(type, _list_to_tensor_space(dims, "dims"))

## dims::AbstractSpace
_gen_dimensions(::Ket, dims::AbstractSpace) = Dimensions(dims, Space(1))
_gen_dimensions(::Bra, dims::AbstractSpace) = Dimensions(Space(1), dims)
_gen_dimensions(::Operator, dims::AbstractSpace) = Dimensions(dims, dims) # is endomorphic
_gen_dimensions(::OperatorKet, dims::AbstractSpace) = Dimensions(dims, Space(1))
_gen_dimensions(::OperatorBra, dims::AbstractSpace) = Dimensions(Space(1), dims)
_gen_dimensions(::SuperOperator, dims::AbstractSpace) = Dimensions(dims, dims) # is endomorphic

## dims::DimsListType{T1, T2} : general nested array
_gen_dimensions(::QuantumObjectType, dims::DimsListType{T1, T2}) where {T1, T2} = Dimensions(dims)

# other cases
_gen_dimensions(::QuantumObjectType, dims::Dimensions) = dims
_gen_dimensions(type::QuantumObjectType, dims) = throw(ArgumentError("The argument `dims` with value $dims is not valid for object type $type."))

# functions for getting Float or Complex element type
_float_type(A::AbstractQuantumObject) = _float_type(eltype(A))
_complex_float_type(A::AbstractQuantumObject) = _complex_float_type(eltype(A))
