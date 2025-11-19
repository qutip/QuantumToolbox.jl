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

function check_dimensions(dimensions_list::NTuple{N, AbstractDimensions}) where {N}
    allequal(dimensions_list) ||
        throw(DimensionMismatch("The quantum objects should have the same Hilbert `dimensions`."))
    return nothing
end
check_dimensions(Qobj_tuple::NTuple{N, AbstractQuantumObject}) where {N} =
    check_dimensions(getfield.(Qobj_tuple, :dimensions))
check_dimensions(A::AbstractQuantumObject...) = check_dimensions(A)

_check_QuantumObject(
    type::ObjType,
    dimensions::GeneralProductDimensions,
    m::Int,
    n::Int,
) where {ObjType <: Union{Ket, Bra, SuperOperator, OperatorBra, OperatorKet}} = throw(
    DomainError(
        _get_dims_string(dimensions),
        "The given `dims` is not compatible with type = $type, should be a single list of integers.",
    ),
)

function _check_QuantumObject(type::Ket, dimensions::ProductDimensions, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Ket"))
    (prod(dimensions) != m) && throw(
        DimensionMismatch("Ket with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n))."),
    )
    return nothing
end

function _check_QuantumObject(type::Bra, dimensions::ProductDimensions, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Bra"))
    (prod(dimensions) != n) && throw(
        DimensionMismatch("Bra with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n))."),
    )
    return nothing
end

function _check_QuantumObject(type::Operator, dimensions::ProductDimensions, m::Int, n::Int)
    L = prod(dimensions)
    (L == m == n) || throw(
        DimensionMismatch(
            "Operator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::Operator, dimensions::GeneralProductDimensions, m::Int, n::Int)
    ((m == 1) || (n == 1)) && throw(DomainError((m, n), "The size of the array is not compatible with Operator"))
    ((prod(dimensions.to) != m) || (prod(dimensions.from) != n)) && throw(
        DimensionMismatch(
            "Operator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::SuperOperator, dimensions::ProductDimensions, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with SuperOperator"))
    (prod(dimensions) != sqrt(m)) && throw(
        DimensionMismatch(
            "SuperOperator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorKet, dimensions::ProductDimensions, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorKet"))
    (prod(dimensions) != sqrt(m)) && throw(
        DimensionMismatch(
            "OperatorKet with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorBra, dimensions::ProductDimensions, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorBra"))
    (prod(dimensions) != sqrt(n)) && throw(
        DimensionMismatch(
            "OperatorBra with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
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

# this returns `to` in GeneralProductDimensions representation
get_dimensions_to(A::AbstractQuantumObject{Ket, <:ProductDimensions}) = A.dimensions.to
get_dimensions_to(A::AbstractQuantumObject{Bra, <:ProductDimensions}) = space_one_list(A.dimensions.to)
get_dimensions_to(A::AbstractQuantumObject{Operator, <:ProductDimensions}) = A.dimensions.to
get_dimensions_to(A::AbstractQuantumObject{Operator, <:GeneralProductDimensions}) = A.dimensions.to
get_dimensions_to(
    A::AbstractQuantumObject{ObjType, <:ProductDimensions},
) where {ObjType <: Union{SuperOperator, OperatorBra, OperatorKet}} = A.dimensions.to

# this returns `from` in GeneralProductDimensions representation
get_dimensions_from(A::AbstractQuantumObject{Ket, <:ProductDimensions}) = space_one_list(A.dimensions.to)
get_dimensions_from(A::AbstractQuantumObject{Bra, <:ProductDimensions}) = A.dimensions.to
get_dimensions_from(A::AbstractQuantumObject{Operator, <:ProductDimensions}) = A.dimensions.to
get_dimensions_from(A::AbstractQuantumObject{Operator, <:GeneralProductDimensions}) = A.dimensions.from
get_dimensions_from(
    A::AbstractQuantumObject{ObjType, <:ProductDimensions},
) where {ObjType <: Union{SuperOperator, OperatorBra, OperatorKet}} = A.dimensions.to

# this creates a list of HilbertSpace(1), it is used to generate `from` for Ket, and `to` for Bra
space_one_list(dimensions::NTuple{N, AbstractSpace}) where {N} =
    ntuple(i -> HilbertSpace(1), Val(sum(length, dimensions)))

# functions for getting Float or Complex element type
_float_type(A::AbstractQuantumObject) = _float_type(eltype(A))
_complex_float_type(A::AbstractQuantumObject) = _complex_float_type(eltype(A))
