export AbstractDimensions, Dimensions#, CompoundDimensions

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

# this creates a list of Space(1), it's used to generate `from` for Ket, and `` for Bra)
# oneDimensions(N::Int) = Dimensions(SVector{N,AbstractSpace}(ntuple(i -> Space(1), Val(N))))

_gen_dims(dims::Union{AbstractVector{T},NTuple{N,T}}) where {T<:Integer,N} = Dimensions(dims)
_gen_dims(dims::Int) = Dimensions(SVector{1,AbstractSpace}(Space(dims)))
_gen_dims(dims::AbstractDimensions) = dims
_gen_dims(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

#= struct CompoundDimensions{N} <: AbstractDimensions{N}
    # note that the number `N` should be the same for both `to` and `from`
    to::SVector{N,AbstractSpace}   # space acting on the left
    from::SVector{N,AbstractSpace} # space acting on the right
end
Base.show(io::IO, D::CompoundDimensions) = print(io, "[", D.to, ", ", D.from, "]")

function CompoundDimensions(to::Union{AbstractVector{T},NTuple{N1,T}}, from::Union{AbstractVector{T},NTuple{N2,T}}) where {T<:Integer,N1,N2}
    _non_static_array_warning("dims", to)
    _non_static_array_warning("dims", from)

    L1 = length(to)
    L2 = length(from)
    ((L1 > 0) && (L1 == L2)) || throw(DomainError((to, from), "The arguments `to` and `from` must be in the same length and have at least one element."))

    return CompoundDimensions(SVector{L1,AbstractSpace}(Space.(to)), SVector{L1,AbstractSpace}(Space.(from)))
end
CompoundDimensions(to::Int, from::Int) = CompoundDimensions(SVector{1,Int}(to), SVector{1,Int}(from))
CompoundDimensions(::KetQuantumObject, dims::Dimensions) = CompoundDimensions(dims, oneDimensions(length(dims)))
CompoundDimensions(::BraQuantumObject, dims::Dimensions) = CompoundDimensions(oneDimensions(length(dims)), dims)
CompoundDimensions(::OperatorQuantumObject, dims::Dimensions) = CompoundDimensions(dims, dims)
CompoundDimensions(::OperatorQuantumObject, dims::CompoundDimensions) = dims =#

Base.length(dims::AbstractDimensions) = length(dims.to)

Base.prod(dims::Dimensions) = prod(dims.to)
Base.prod(spaces::SVector{1,AbstractSpace}) = spaces[1].size # for `Dimensions.to` has only a single Space

LinearAlgebra.transpose(dims::Dimensions) = dims
# LinearAlgebra.transpose(dims::CompoundDimensions) = CompoundDimensions(dims.from, dims.to) # switch `to` and `from`

LinearAlgebra.kron(Adims::Dimensions, Bdims::Dimensions) = Dimensions(vcat(Adims.to, Bdims.to))
# LinearAlgebra.kron(Adims::CompoundDimensions, Bdims::CompoundDimensions) = CompoundDimensions(vcat(Adims.to, Bdims.to), vcat(Adims.from, Bdims.from))