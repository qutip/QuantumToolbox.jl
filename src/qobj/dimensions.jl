#=
This file defines the Dimensions structures, which can describe composite Hilbert spaces.
=#

export AbstractDimensions, Dimensions, GeneralDimensions

abstract type AbstractDimensions{M,N} end

@doc raw"""
    struct ProductDimensions{M,N,T1<:Tuple,T2<:Union{<:Tuple, Nothing}} <: AbstractDimensions{M,N}
        to::T1
        from::T2
    end

A structure that describes the left-hand side (`to`) and right-hand side (`from`) dimensions of a quantum object.
"""
struct ProductDimensions{M,N,T1<:Tuple,T2<:Union{<:Tuple,Nothing}} <: AbstractDimensions{M,N}
    to::T1   # space acting on the left
    from::T2 # space acting on the right

    # make sure the elements in the tuple are all AbstractSpace
    GeneralDimensions(to::NTuple{M,AbstractSpace}, from::Union{NTuple{N,AbstractSpace},Nothing}) where {M,N} =
        new{M,N,typeof(to),typeof(from)}(to, from)
end
function ProductDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    to = _dims_tuple_of_space(dims[1])
    from = isnothing(dims[2]) ? nothing : _dims_tuple_of_space(dims[2])

    return ProductDimensions(to, from)
end

ProductDimensions(dims::Union{Int,AbstractSpace}) = ProductDimensions((dims,), nothing)
ProductDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = ProductDimensions(dims, nothing)
ProductDimensions(to::NTuple{M,Int}, from::Union{NTuple{N,Int},Nothing}) where {M,N} = ProductDimensions((to, from))
ProductDimensions(dims::ProductDimensions) = dims

# obtain dims in the type of SVector with integers
dimensions_to_dims(dimensions::NTuple{N,AbstractSpace}) where {N} = SVector{N}(map(dimensions_to_dims, dimensions))
dimensions_to_dims(dimensions::ProductDimensions) =
    SVector{2}(dimensions_to_dims(dimensions.to), dimensions_to_dims(dimensions.from))
dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

hilbert_dimensions_to_size(dimensions::ProductDimensions) =
    (hilbert_dimensions_to_size(dimensions.to), hilbert_dimensions_to_size(dimensions.from))
hilbert_dimensions_to_size(dimensions::NTuple{N,AbstractSpace}) where {N} = prod(hilbert_dimensions_to_size, dimensions)
hilbert_dimensions_to_size(::Nothing) = nothing

liouvillian_dimensions_to_size(dimensions::ProductDimensions) =
    (liouvillian_dimensions_to_size(dimensions.to), liouvillian_dimensions_to_size(dimensions.from))
liouvillian_dimensions_to_size(dimensions::NTuple{N,AbstractSpace}) where {N} =
    prod(liouvillian_dimensions_to_size, dimensions)
liouvillian_dimensions_to_size(::Nothing) = nothing

Base.length(::AbstractDimensions{N}) where {N} = N

Base.transpose(dimensions::ProductDimensions) = ProductDimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::AbstractDimensions) = transpose(dimensions)

Base.:(==)(dim1::ProductDimensions, dim2::ProductDimensions) = (dim1.to == dim2.to) && (dim1.from == dim2.from)

function _dims_tuple_of_space(dims::NTuple{N,AbstractSpace}) where {N}
    _non_static_array_warning("dims", dims)

    N > 0 || throw(DomainError(N, "The length of `dims` must be larger or equal to 1."))

    return ntuple(dim -> HilbertSpace(dims[dim]), Val(N))
end

# this is used to show `dims` for Qobj and QobjEvo
function _get_dims_string(dimensions::ProductDimensions)
    dims = dimensions_to_dims(dimensions)
    isnothing(dims[2]) && return string(dims[1])
    return "[$(string(dims[1])), $(string(dims[2]))]"
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing
