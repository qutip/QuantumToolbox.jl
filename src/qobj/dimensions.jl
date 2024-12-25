export AbstractDimensions, Dimensions, CompoundDimensions

abstract type AbstractDimensions{N} end

struct Dimensions{N} <: AbstractDimensions{N}
    to::SVector{N,<:AbstractSpace}
end
function Dimensions(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    return Dimensions{L}(SVector{L,Space}(Space.(dims)))
end
Dimensions(dims::Int) = Dimensions(SVector{1,Int}(dims))
Dimensions(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

Base.show(io::IO, D::Dimensions) = print(io, D.to)

# this creates a list of Space(1), it's used to generate `from` for Ket, and `to` for Bra)
oneDimensions(N::Int) = Dimensions(SVector{N,Space}(ntuple(i -> Space(1), Val(N))))

struct CompoundDimensions{N} <: AbstractDimensions{N}
    # note that the number `N` should be the same for both `to` and `from`
    to::SVector{N,<:AbstractSpace}   # space acting on the left
    from::SVector{N,<:AbstractSpace} # space acting on the right
end
function CompoundDimensions(to::Union{AbstractVector{T},NTuple{N1,T}}, from::Union{AbstractVector{T},NTuple{N2,T}}) where {T<:Integer,N1,N2}
    _non_static_array_warning("dims", to)
    _non_static_array_warning("dims", from)

    L1 = length(to)
    L2 = length(from)
    ((L1 > 0) && (L1 == L2)) || throw(DomainError((L1, L2), "The length of the arguments `to` and `from` must be in the same length and have at least one element."))

    return CompoundDimensions{L1}(SVector{L1,Space}(Space.(to)), SVector{L1,Space}(Space.(from)))
end
CompoundDimensions(to::Int, from::Int) = CompoundDimensions(SVector{1,Int}(to), SVector{1,Int}(from))
CompoundDimensions(::KetQuantumObject, dims::Dimensions{N}) where {N} = CompoundDimensions{N}(dims, oneDimensions(length(dims)))
CompoundDimensions(::BraQuantumObject, dims::Dimensions{N}) where {N} = CompoundDimensions{N}(oneDimensions(length(dims)), dims)
CompoundDimensions(::OperatorQuantumObject, dims::Dimensions{N}) where {N} = CompoundDimensions{N}(dims, dims)
CompoundDimensions(::OperatorQuantumObject, dims::CompoundDimensions) = dims

Base.show(io::IO, D::CompoundDimensions) = print(io, "[", D.to, ", ", D.from, "]")

_gen_dims(dims::AbstractDimensions) = dims
_gen_dims(dims::Any) = Dimensions(dims)

Base.length(dims::AbstractDimensions) = length(dims.to)

Base.prod(dims::Dimensions) = prod(dims.to)
Base.prod(spaces::SVector{1,<:AbstractSpace}) = spaces[1].size # for `Dimensions.to` has only a single Space

LinearAlgebra.transpose(dims::Dimensions) = dims
LinearAlgebra.transpose(dims::CompoundDimensions) = CompoundDimensions(dims.from, dims.to) # switch `to` and `from`

LinearAlgebra.kron(Adims::Dimensions{NA}, Bdims::Dimensions{NB}) where {NA,NB} = Dimensions{NA+NB}(vcat(Adims.to, Bdims.to))
LinearAlgebra.kron(Adims::CompoundDimensions{NA}, Bdims::CompoundDimensions{NB}) where {NA,NB} = CompoundDimensions{NA+NB}(vcat(Adims.to, Bdims.to), vcat(Adims.from, Bdims.from))