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
    abstract type AbstractQuantumObject{DataType,ObjType,DimType}

Abstract type for all quantum objects like [`QuantumObject`](@ref) and [`QuantumObjectEvolution`](@ref).

# Example
```jldoctest
julia> sigmax() isa AbstractQuantumObject
true
```
"""
abstract type AbstractQuantumObject{DataType,ObjType,DimType} end

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
    shape(A::AbstractQuantumObject)
    shape(A::AbstractQuantumObject, idx::Int)

Returns a tuple containing each dimensions of the array in the [`AbstractQuantumObject`](@ref).

Optionally, you can specify an index (`idx`) to just get the corresponding dimension of the array.

!!! note
    `shape` is a synonym of `size`.
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
    isequal(A.type, B.type) && isequal(A.dimensions, B.dimensions) && isequal(A.data, B.data)
Base.isapprox(A::AbstractQuantumObject, B::AbstractQuantumObject; kwargs...) =
    isequal(A.type, B.type) && isequal(A.dimensions, B.dimensions) && isapprox(A.data, B.data; kwargs...)
Base.:(==)(A::AbstractQuantumObject, B::AbstractQuantumObject) =
    (A.type == B.type) && (A.dimensions == B.dimensions) && (A.data == B.data)

function check_dimensions(dimensions_list::NTuple{N,AbstractDimensions}) where {N}
    allequal(dimensions_list) ||
        throw(DimensionMismatch("The quantum objects should have the same Hilbert `dimensions`."))
    return nothing
end
check_dimensions(Qobj_tuple::NTuple{N,AbstractQuantumObject}) where {N} =
    check_dimensions(getfield.(Qobj_tuple, :dimensions))
check_dimensions(A::AbstractQuantumObject...) = check_dimensions(A)

_check_QuantumObject(
    type::ObjType,
    dimensions::GeneralDimensions,
    m::Int,
    n::Int,
) where {
    ObjType<:Union{
        KetQuantumObject,
        BraQuantumObject,
        SuperOperatorQuantumObject,
        OperatorBraQuantumObject,
        OperatorKetQuantumObject,
    },
} = throw(
    DomainError(
        _get_dims_string(dimensions),
        "The given `dims` is not compatible with type = $type, should be an `Operator`.",
    ),
)

function _check_QuantumObject(type::KetQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Ket"))
    (prod(dimensions) != m) && throw(
        DimensionMismatch("Ket with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n))."),
    )
    return nothing
end

function _check_QuantumObject(type::BraQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with Bra"))
    (prod(dimensions) != n) && throw(
        DimensionMismatch("Bra with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n))."),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    L = prod(dimensions)
    (L == m == n) || throw(
        DimensionMismatch(
            "Operator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorQuantumObject, dimensions::GeneralDimensions, m::Int, n::Int)
    ((m == 1) || (n == 1)) && throw(DomainError((m, n), "The size of the array is not compatible with Operator"))
    ((prod(dimensions.to) != m) || (prod(dimensions.from) != n)) && throw(
        DimensionMismatch(
            "Operator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::SuperOperatorQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    (m != n) && throw(DomainError((m, n), "The size of the array is not compatible with SuperOperator"))
    (prod(dimensions) != sqrt(m)) && throw(
        DimensionMismatch(
            "SuperOperator with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorKetQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    (n != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorKet"))
    (prod(dimensions) != sqrt(m)) && throw(
        DimensionMismatch(
            "OperatorKet with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

function _check_QuantumObject(type::OperatorBraQuantumObject, dimensions::Dimensions, m::Int, n::Int)
    (m != 1) && throw(DomainError((m, n), "The size of the array is not compatible with OperatorBra"))
    (prod(dimensions) != sqrt(n)) && throw(
        DimensionMismatch(
            "OperatorBra with dims = $(_get_dims_string(dimensions)) does not fit the array size = $((m, n)).",
        ),
    )
    return nothing
end

Base.getproperty(A::AbstractQuantumObject, key::Symbol) = getproperty(A, makeVal(key))

# support `AbstractQuantumObject.dims`
Base.getproperty(A::AbstractQuantumObject, ::Val{:dims}) = dimensions_to_dims(getfield(A, :dimensions))
Base.getproperty(A::AbstractQuantumObject, ::Val{K}) where {K} = getfield(A, K)

# this returns `to` in GeneralDimensions representation
get_dimensions_to(A::AbstractQuantumObject{DT,KetQuantumObject,Dimensions{N}}) where {DT,N} = A.dimensions.to
get_dimensions_to(A::AbstractQuantumObject{DT,BraQuantumObject,Dimensions{N}}) where {DT,N} = space_one_list(N)
get_dimensions_to(A::AbstractQuantumObject{DT,OperatorQuantumObject,Dimensions{N}}) where {DT,N} = A.dimensions.to
get_dimensions_to(A::AbstractQuantumObject{DT,OperatorQuantumObject,GeneralDimensions{N}}) where {DT,N} =
    A.dimensions.to
get_dimensions_to(
    A::AbstractQuantumObject{DT,ObjType,Dimensions{N}},
) where {DT,ObjType<:Union{SuperOperatorQuantumObject,OperatorBraQuantumObject,OperatorKetQuantumObject},N} =
    A.dimensions.to

# this returns `from` in GeneralDimensions representation
get_dimensions_from(A::AbstractQuantumObject{DT,KetQuantumObject,Dimensions{N}}) where {DT,N} = space_one_list(N)
get_dimensions_from(A::AbstractQuantumObject{DT,BraQuantumObject,Dimensions{N}}) where {DT,N} = A.dimensions.to
get_dimensions_from(A::AbstractQuantumObject{DT,OperatorQuantumObject,Dimensions{N}}) where {DT,N} = A.dimensions.to
get_dimensions_from(A::AbstractQuantumObject{DT,OperatorQuantumObject,GeneralDimensions{N}}) where {DT,N} =
    A.dimensions.from
get_dimensions_from(
    A::AbstractQuantumObject{DT,ObjType,Dimensions{N}},
) where {DT,ObjType<:Union{SuperOperatorQuantumObject,OperatorBraQuantumObject,OperatorKetQuantumObject},N} =
    A.dimensions.to

# functions for getting Float or Complex element type
_FType(A::AbstractQuantumObject) = _FType(eltype(A))
_CType(A::AbstractQuantumObject) = _CType(eltype(A))
