#=
This file defines the Dimensions structures, which can describe composite Hilbert spaces.
=#

export AbstractDimensions, Dimensions, GeneralDimensions

abstract type AbstractDimensions{M,N} end

@doc raw"""
    struct ProductDimensions{M,N,T1<:Tuple,T2<:Tuple} <: AbstractDimensions{M,N}
        to::T1
        from::T2
    end

A structure that describes the left-hand side (`to`) and right-hand side (`from`) dimensions of a quantum object.
"""
struct ProductDimensions{M,N,T1<:Tuple,T2<:Tuple} <: AbstractDimensions{M,N}
    to::T1   # space acting on the left
    from::T2 # space acting on the right

    # make sure the elements in the tuple are all AbstractSpace
    GeneralDimensions(to::NTuple{M,AbstractSpace}, from::NTuple{N,AbstractSpace}) where {M,N} =
        new{M,N,typeof(to),typeof(from)}(to, from)
end
function ProductDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    _non_static_array_warning("dims[1]", dims[1])
    _non_static_array_warning("dims[2]", dims[2])

    L1 = length(dims[1])
    L2 = length(dims[2])
    (L1 > 0) || throw(DomainError(L1, "The length of `dims[1]` must be larger or equal to 1."))
    (L2 > 0) || throw(DomainError(L2, "The length of `dims[2]` must be larger or equal to 1."))

    to = ntuple(i -> HilbertSpace(dims[1][i]), Val(L1))
    from = ntuple(i -> HilbertSpace(dims[2][i]), Val(L2))

    return ProductDimensions(to, from)
end

ProductDimensions(dims::Union{Int, AbstractSpace}) = ProductDimensions((dims,), (dims,))
ProductDimensions(to::NTuple{M,Int}, from::NTuple{N,Int}) where {M,N} = ProductDimensions((to, from))
ProductDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = ProductDimensions(dims, dims)
ProductDimensions(dims::ProductDimensions) = dims

# obtain dims in the type of SVector with integers
dimensions_to_dims(dimensions::NTuple{N,AbstractSpace}) where {N} = SVector{N}(map(dimensions_to_dims, dimensions))
dimensions_to_dims(dimensions::ProductDimensions) =
    SVector{2}(dimensions_to_dims(dimensions.to), dimensions_to_dims(dimensions.from))
dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

hilbert_dimensions_to_size(dimensions::ProductDimensions) = (prod(hilbert_dimensions_to_size, dimensions.to), prod(hilbert_dimensions_to_size, dimensions.from))
liouvillian_dimensions_to_size(dimensions::ProductDimensions) =
    (prod(liouvillian_dimensions_to_size, dimensions.to), prod(liouvillian_dimensions_to_size, dimensions.from))


Base.length(::AbstractDimensions{N}) where {N} = N

Base.transpose(dimensions::ProductDimensions) = ProductDimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::AbstractDimensions) = transpose(dimensions)

# this is used to show `dims` for Qobj and QobjEvo
_get_dims_string(dimensions::Dimensions) = string(dimensions_to_dims(dimensions))
function _get_dims_string(dimensions::ProductDimensions)
    dims = dimensions_to_dims(dimensions)
    dims[1] == dims[2] && return string(dims[1])
    return "[$(string(dims[1])), $(string(dims[2]))]"
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing

Base.:(==)(dim1::ProductDimensions, dim2::ProductDimensions) = (dim1.to == dim2.to) && (dim1.from == dim2.from)
