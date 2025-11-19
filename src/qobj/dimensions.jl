#=
This file defines the Dimensions structures, which can describe composite Hilbert spaces.
=#

export AbstractDimensions, ProductDimensions

abstract type AbstractDimensions{M,N} end

@doc raw"""
    struct ProductDimensions{M,N,T1<:Tuple,T2<:Union{<:Tuple, Nothing}} <: AbstractDimensions{M,N}
        to::T1
        from::T2
    end

A structure that describes the left-hand side (`to`) and right-hand side (`from`) dimensions of a quantum object.

The `from` field can be different from `nothing` only for non-square [`Operator`](@ref) and [`SuperOperator`](@ref) quantum objects.
"""
struct ProductDimensions{M,N,T1<:Tuple,T2<:Union{<:Tuple,Nothing}} <: AbstractDimensions{M,N}
    to::T1   # space acting on the left
    from::T2 # space acting on the right

    function ProductDimensions(to::Union{AbstractVector,NTuple}, from::Union{NTuple,Nothing})
        M = length(to)
        N = isnothing(from) ? M : length(from)

        _non_static_array_warning("dims", to)
        isnothing(from) || _non_static_array_warning("dims", from)

        to_space = _dims_tuple_of_space(to)
        from_space = _dims_tuple_of_space(from)

        new{M,N,typeof(to_space),typeof(from_space)}(to_space, from_space)
    end
end
function ProductDimensions(dims::Union{AbstractVector,Tuple})
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    return ProductDimensions(dims[1], dims[2])
end

ProductDimensions(dims::Union{Int,AbstractSpace}) = ProductDimensions((dims,), nothing)
ProductDimensions(dims::Union{AbstractVector{<:Integer},NTuple{N,Integer}}) where {N} = ProductDimensions(dims, nothing)
ProductDimensions(dims::Union{AbstractVector{<:AbstractSpace},NTuple{N,AbstractSpace}}) where {N} = ProductDimensions(dims, nothing)
ProductDimensions(dims::ProductDimensions) = dims

# obtain dims in the type of SVector with integers
dimensions_to_dims(dimensions::NTuple{N,AbstractSpace}) where {N} = vcat(map(dimensions_to_dims, dimensions)...)
function dimensions_to_dims(dimensions::ProductDimensions)
    dims_to = dimensions_to_dims(dimensions.to)
    isnothing(dimensions.from) && return dims_to
    dims_from = dimensions_to_dims(dimensions.from)
    return SVector{2}(dims_to, dims_from)
end
dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

hilbert_dimensions_to_size(dimensions::ProductDimensions) =
    (hilbert_dimensions_to_size(dimensions.to), hilbert_dimensions_to_size(dimensions.from))
hilbert_dimensions_to_size(dimensions::NTuple{N,AbstractSpace}) where {N} = prod(hilbert_dimensions_to_size, dimensions)
hilbert_dimensions_to_size(dim::Int) = dim
hilbert_dimensions_to_size(dimensions::Union{AbstractVector, NTuple{N,Integer}}) where {N} = prod(dimensions)
hilbert_dimensions_to_size(::Nothing) = nothing

liouvillian_dimensions_to_size(dimensions::ProductDimensions) =
    (liouvillian_dimensions_to_size(dimensions.to), liouvillian_dimensions_to_size(dimensions.from))
liouvillian_dimensions_to_size(dimensions::NTuple{N,AbstractSpace}) where {N} =
    prod(liouvillian_dimensions_to_size, dimensions)
liouvillian_dimensions_to_size(::Nothing) = nothing

Base.length(::AbstractDimensions{N}) where {N} = N

Base.transpose(dimensions::ProductDimensions) = ProductDimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.transpose(dimensions::ProductDimensions{M,N,T1,Nothing}) where {M,N,T1<:Tuple} = dimensions
Base.adjoint(dimensions::AbstractDimensions) = transpose(dimensions)

Base.:(==)(dim1::ProductDimensions, dim2::ProductDimensions) = (dim1.to == dim2.to) && (dim1.from == dim2.from)

_dims_tuple_of_space(dims::NTuple{N,AbstractSpace}) where {N} = dims
function _dims_tuple_of_space(dims::Union{AbstractVector{<:Integer},SVector{M, <:Integer}, NTuple{M, Integer}}) where {M}
    _non_static_array_warning("dims", dims)

    N = length(dims)
    N > 0 || throw(DomainError(N, "The length of `dims` must be larger or equal to 1."))

    return ntuple(dim -> HilbertSpace(dims[dim]), Val(N))
end
_dims_tuple_of_space(::Nothing) = nothing

_gen_dimensions(dims::AbstractDimensions) = dims
_gen_dimensions(dims) = ProductDimensions(dims)

# this is used to show `dims` for Qobj and QobjEvo
function _get_dims_string(dimensions::ProductDimensions)
    dims = dimensions_to_dims(dimensions)
    isnothing(dimensions.from) && return string(dims)
    return "[$(string(dims[1])), $(string(dims[2]))]"
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing

space_one_list(dimensions::NTuple{N,AbstractSpace}) where {N} =
    ntuple(i -> one(dimensions[i]), Val(sum(length, dimensions)))
