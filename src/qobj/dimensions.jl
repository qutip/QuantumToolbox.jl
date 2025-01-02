export AbstractDimensions, Dimensions, CompoundDimensions

abstract type AbstractDimensions{N} end

# this show function is for printing AbstractDimensions
function Base.show(io::IO, svec::SVector{N,AbstractSpace}) where {N}
    print(io, "[")
    join(io, string.(svec), ", ")
    return print(io, "]")
end

struct Dimensions{N} <: AbstractDimensions{N}
    to::SVector{N,AbstractSpace}
end
function Dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    return Dimensions{L}(SVector{L,AbstractSpace}(Space.(dims)))
end
Dimensions(dims::Int) = Dimensions(SVector{1,Int}(dims))
Dimensions(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

Base.show(io::IO, D::Dimensions) = print(io, D.to)

struct CompoundDimensions{N} <: AbstractDimensions{N}
    # note that the number `N` should be the same for both `to` and `from`
    to::SVector{N,AbstractSpace}   # space acting on the left
    from::SVector{N,AbstractSpace} # space acting on the right
end
function CompoundDimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    _non_static_array_warning("dims[1]", dims[1])
    _non_static_array_warning("dims[2]", dims[2])

    L1 = length(dims[1])
    L2 = length(dims[2])
    ((L1 > 0) && (L1 == L2)) || throw(
        DomainError(
            (L1, L2),
            "The length of the arguments `dims[1]` and `dims[2]` must be in the same length and have at least one element.",
        ),
    )

    return CompoundDimensions{L1}(
        SVector{L1,AbstractSpace}(Space.(dims[1])),
        SVector{L1,AbstractSpace}(Space.(dims[2])),
    )
end

Base.show(io::IO, D::CompoundDimensions) = print(io, "[", D.to, ", ", D.from, "]")

_gen_dims(dims::AbstractDimensions) = dims
_gen_dims(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = Dimensions(dims)
_gen_dims(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Union{AbstractVector,NTuple},N} =
    CompoundDimensions(dims)
_gen_dims(dims::Any) = Dimensions(dims)

# obtain dims in the type of SVector with integers
dims_to_list(dimsvec::SVector{N,AbstractSpace}) where {N} = SVector{N,Int}(ntuple(i -> dimsvec[i].size, Val(N)))
dims_to_list(dims::Dimensions) = dims_to_list(dims.to)
dims_to_list(dims::CompoundDimensions) = SVector{2}(dims_to_list(dims.to), dims_to_list(dims.from))

Base.:(==)(vect::AbstractVector{T}, dims::AbstractDimensions) where {T} = vect == dims_to_list(dims)
Base.:(==)(dims::AbstractDimensions, vect::AbstractVector{T}) where {T} = vect == dims

Base.length(::AbstractDimensions{N}) where {N} = N

# need to specify return type `Int` for `_get_space_size`
# otherwise the type of `prod(::Dimensions)` will be unstable
_get_space_size(s::AbstractSpace)::Int = s.size
Base.prod(dims::Dimensions) = prod(dims.to)
Base.prod(spaces::SVector{N,AbstractSpace}) where {N} = prod(_get_space_size, spaces)

LinearAlgebra.transpose(dims::Dimensions) = dims
LinearAlgebra.transpose(dims::CompoundDimensions) = CompoundDimensions(dims.from, dims.to) # switch `to` and `from`
