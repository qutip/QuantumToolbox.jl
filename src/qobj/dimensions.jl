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
function CompoundDimensions(
    to::Union{AbstractVector{T},NTuple{N1,T}},
    from::Union{AbstractVector{T},NTuple{N2,T}},
) where {T<:Integer,N1,N2}
    _non_static_array_warning("dims", to)
    _non_static_array_warning("dims", from)

    L1 = length(to)
    L2 = length(from)
    ((L1 > 0) && (L1 == L2)) || throw(
        DomainError(
            (L1, L2),
            "The length of the arguments `to` and `from` must be in the same length and have at least one element.",
        ),
    )

    return CompoundDimensions{L1}(SVector{L1,AbstractSpace}(Space.(to)), SVector{L1,AbstractSpace}(Space.(from)))
end
CompoundDimensions(to::Int, from::Int) = CompoundDimensions(SVector{1,Int}(to), SVector{1,Int}(from))

Base.show(io::IO, D::CompoundDimensions) = print(io, "[", D.to, ", ", D.from, "]")

_gen_dims(dims::AbstractDimensions) = dims
_gen_dims(dims::Any) = Dimensions(dims)

# obtain dims in the type of SVector with integers
dims_to_list(dimsvec::SVector{N,AbstractSpace}) where {N} = SVector{N,Int}(ntuple(i -> dimsvec[i].size, Val(N)))
dims_to_list(dims::Dimensions) = dims_to_list(dims.to)
dims_to_list(dims::CompoundDimensions) = SVector{2}(dims_to_list(dims.to), dims_to_list(dims.from))

Base.:(==)(vect::AbstractVector{T}, dims::AbstractDimensions) where {T} = vect == dims_to_list(dims)
Base.:(==)(dims::AbstractDimensions, vect::AbstractVector{T}) where {T} = vect == dims

Base.length(::AbstractDimensions{N}) where {N} = N

Base.prod(dims::Dimensions) = prod(dims.to)

LinearAlgebra.transpose(dims::Dimensions) = dims
LinearAlgebra.transpose(dims::CompoundDimensions) = CompoundDimensions(dims.from, dims.to) # switch `to` and `from`

LinearAlgebra.kron(Adims::Dimensions{NA}, Bdims::Dimensions{NB}) where {NA,NB} =
    Dimensions{NA + NB}(vcat(Adims.to, Bdims.to))
LinearAlgebra.kron(Adims::CompoundDimensions{NA}, Bdims::CompoundDimensions{NB}) where {NA,NB} =
    CompoundDimensions{NA + NB}(vcat(Adims.to, Bdims.to), vcat(Adims.from, Bdims.from))
