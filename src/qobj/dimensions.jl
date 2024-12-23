export AbstractDimensions, Dimensions

abstract type AbstractDimensions{N} end

struct Dimensions{N} <: AbstractDimensions{N}
    to::SVector{N,AbstractSpace}
end
Base.show(io::IO, D::Dimensions) = print(io, D.to)

function Dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    return Dimensions(SVector{L,AbstractSpace}(Space.(dims)))
end

_gen_dims(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = Dimensions(dims)
_gen_dims(dims::Int) = Dimensions(SVector{1,AbstractSpace}(Space(dims)))
_gen_dims(dims::AbstractDimensions) = dims
_gen_dims(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

Base.prod(dims::Dimensions) = prod(dims.to)
Base.prod(spaces::SVector{1,AbstractSpace}) = spaces[1].size # for `Dimensions.to` has only a single Space

LinearAlgebra.transpose(dims::Dimensions) = dims

struct CompoundDimensions{N} <: AbstractDimensions{N}
    # note that the number `N` should be the same for both `to` and `from`
    to::SVector{N,AbstractSpace}   # space acting on the left
    from::SVector{N,AbstractSpace} # space acting on the right
end
Base.show(io::IO, D::CompoundDimensions) = print(io, "[", D.to, ", ", D.from, "]")

LinearAlgebra.transpose(dims::CompoundDimensions) = CompoundDimensions(dims.from, dims.to) # switch `to` and `from`