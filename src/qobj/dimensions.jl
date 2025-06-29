#=
This file defines the Dimensions structures, which can describe composite Hilbert spaces.
=#

export AbstractDimensions, Dimensions, GeneralDimensions

abstract type AbstractDimensions{M,N} end

@doc raw"""
    struct Dimensions{N,T<:Tuple} <: AbstractDimensions{N, N}
        to::T
    end

A structure that describes the Hilbert [`Space`](@ref) of each subsystems.
"""
struct Dimensions{N,T<:Tuple} <: AbstractDimensions{N,N}
    to::T

    # make sure the elements in the tuple are all AbstractSpace
    Dimensions(to::NTuple{N,AbstractSpace}) where {N} = new{N,typeof(to)}(to)
end
function Dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    return Dimensions(Tuple(Space.(dims)))
end
Dimensions(dims::Int) = Dimensions(Space(dims))
Dimensions(dims::DimType) where {DimType<:AbstractSpace} = Dimensions((dims,))
Dimensions(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

@doc raw"""
    struct GeneralDimensions{N,T1<:Tuple,T2<:Tuple} <: AbstractDimensions{N}
        to::T1
        from::T2
    end

A structure that describes the left-hand side (`to`) and right-hand side (`from`) Hilbert [`Space`](@ref) of an [`Operator`](@ref).
"""
struct GeneralDimensions{M,N,T1<:Tuple,T2<:Tuple} <: AbstractDimensions{M,N}
    # note that the number `N` should be the same for both `to` and `from`
    to::T1   # space acting on the left
    from::T2 # space acting on the right

    # make sure the elements in the tuple are all AbstractSpace
    GeneralDimensions(to::NTuple{M,AbstractSpace}, from::NTuple{N,AbstractSpace}) where {M,N} =
        new{M,N,typeof(to),typeof(from)}(to, from)
end
function GeneralDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    _non_static_array_warning("dims[1]", dims[1])
    _non_static_array_warning("dims[2]", dims[2])

    L1 = length(dims[1])
    L2 = length(dims[2])
    (L1 > 0) || throw(DomainError(L1, "The length of `dims[1]` must be larger or equal to 1."))
    (L2 > 0) || throw(DomainError(L2, "The length of `dims[2]` must be larger or equal to 1."))

    return GeneralDimensions(Tuple(Space.(dims[1])), Tuple(Space.(dims[2])))
end

_gen_dimensions(dims::AbstractDimensions) = dims
_gen_dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = Dimensions(dims)
_gen_dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N} =
    GeneralDimensions(dims)
_gen_dimensions(dims::Any) = Dimensions(dims)

# obtain dims in the type of SVector with integers
dimensions_to_dims(dimensions::NTuple{N,AbstractSpace}) where {N} = vcat(map(dimensions_to_dims, dimensions)...)
dimensions_to_dims(dimensions::Dimensions) = dimensions_to_dims(dimensions.to)
dimensions_to_dims(dimensions::GeneralDimensions) =
    SVector{2}(dimensions_to_dims(dimensions.to), dimensions_to_dims(dimensions.from))

dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

Base.length(::AbstractDimensions{N}) where {N} = N

# need to specify return type `Int` for `_get_space_size`
# otherwise the type of `prod(::Dimensions)` will be unstable
_get_space_size(s::AbstractSpace)::Int = s.size
Base.prod(dims::Dimensions) = prod(dims.to)
Base.prod(spaces::NTuple{N,AbstractSpace}) where {N} = prod(_get_space_size, spaces)

Base.transpose(dimensions::Dimensions) = dimensions
Base.transpose(dimensions::GeneralDimensions) = GeneralDimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::AbstractDimensions) = transpose(dimensions)

# this is used to show `dims` for Qobj and QobjEvo
_get_dims_string(dimensions::Dimensions) = string(dimensions_to_dims(dimensions))
function _get_dims_string(dimensions::GeneralDimensions)
    dims = dimensions_to_dims(dimensions)
    return "[$(string(dims[1])), $(string(dims[2]))]"
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing

Base.:(==)(dim1::Dimensions, dim2::Dimensions) = (dim1.to == dim2.to)
Base.:(==)(dim1::GeneralDimensions, dim2::GeneralDimensions) = (dim1.to == dim2.to) && (dim1.from == dim2.from)
Base.:(==)(dim1::Dimensions, dim2::GeneralDimensions) = false
Base.:(==)(dim1::GeneralDimensions, dim2::Dimensions) = false
