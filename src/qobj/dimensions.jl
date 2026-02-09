#=
This file defines the Dimensions structure and also the following space structures:
    - Space
    - TensorSpace
    - LiouvilleSpace
=#

export Dimensions
export AbstractSpace, Space, TensorSpace, LiouvilleSpace

@doc raw"""
    abstract type AbstractSpace

Abstract type for all space structures.
"""
abstract type AbstractSpace end

##########################################
# Dimensions

@doc raw"""
    struct Dimensions{T1 <: AbstractSpace, T2 <: AbstractSpace}
        to::T1
        from::T2
    end

A structure that embodies the left-hand side (`to`) and right-hand side (`from`) [`AbstractSpace`](@ref) of a quantum object.

The fields `to` and `from` are related to the left (row) and right (column) dimensions, respectively.

# Constructors

- `Dimensions(to, from)`: Create with explicit `to` and `from` [`AbstractSpace`](@ref)
- `Dimensions(dims)`: Create square dimensions where `to == from`
- `Dimensions((to_dims, from_dims))`: Create non-square dimensions from a 2-element tuple of integer vectors

# Examples

```jldoctest
julia> Dimensions(3)  # Single 3-dimensional Hilbert space
Dimensions{1, 1, Tuple{Space}, Tuple{Space}}((Space(3),), (Space(3),))

julia> Dimensions((2, 3))  # Composite 2⊗3 Hilbert space (square)
Dimensions{2, 2, Tuple{HilbertSpace, HilbertSpace}, Tuple{HilbertSpace, HilbertSpace}}((HilbertSpace(2), HilbertSpace(3)), (HilbertSpace(2), HilbertSpace(3)))

julia> Dimensions(((2, 3), (4,)))  # Non-square: maps from 4-dim to 2⊗3=6-dim
Dimensions{2, 1, Tuple{HilbertSpace, HilbertSpace}, Tuple{HilbertSpace}}((HilbertSpace(2), HilbertSpace(3)), (HilbertSpace(4),))
```
"""
struct Dimensions{T1 <: AbstractSpace, T2 <: AbstractSpace}
    to::T1   # space acting on the left (rows)
    from::T2 # space acting on the right (columns)
end

# Square dimensions from integer tuple/vector
function Dimensions(dims::Union{AbstractVector{T}, NTuple{N, T}}) where {T <: Integer, N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    spaces = TensorSpace(Tuple(Space.(dims)))
    return Dimensions(spaces, spaces)
end

# Square dimensions from single integer
Dimensions(dims::Int) = Dimensions(Space(dims))

# Square dimensions from single AbstractSpace
Dimensions(dims::AbstractSpace) = Dimensions(dims, dims)

# Non-square dimensions from 2-element tuple of integer vectors/tuples
function Dimensions(dims::Union{AbstractVector{T}, NTuple{2, T}}) where {T <: Union{AbstractVector, NTuple}}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    _non_static_array_warning("dims[1]", dims[1])
    _non_static_array_warning("dims[2]", dims[2])

    L1 = length(dims[1])
    L2 = length(dims[2])
    (L1 > 0) || throw(DomainError(L1, "The length of `dims[1]` must be larger or equal to 1."))
    (L2 > 0) || throw(DomainError(L2, "The length of `dims[2]` must be larger or equal to 1."))

    return Dimensions(TensorSpace(Tuple(Space.(dims[1]))), TensorSpace(Tuple(Space.(dims[2]))))
end

# Error for invalid input
Dimensions(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

function Base.show(io::IO, d::Dimensions)
    print(io, "Dimensions($(d.to), $(d.from))")
    return nothing
end

# Check if dimensions are square (to == from)
isendomorphism(dimensions::Dimensions) = dimensions.to == dimensions.from

# obtain dims in the type of SVector with integers
function dimensions_to_dims(dimensions::Dimensions)
    to_dims = dimensions_to_dims(dimensions.to)
    from_dims = dimensions_to_dims(dimensions.from)
    return (to_dims, from_dims)
end

dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

"""
    get_size(dimensions::Dimensions)

Returns the matrix dimensions `(m, n)` of a given [`Dimensions`](@ref).

Returns `(m, n)` where `m` is the product of the `dimensions.to`, and `n` is the product of the `dimensions.from`.

If `dimensions` is an `Integer` or a vector/tuple of `Integer`s, it is automatically treated as `Dimensions(dimensions, dimensions)`.
"""
get_size(dimensions::Dimensions) = (get_size(dimensions.to), get_size(dimensions.from))
get_size(dimensions::Union{<:Integer, AbstractVector{<:Integer}, NTuple{N, Integer}}) where {N} = get_size(Dimensions(dimensions))

Base.transpose(dimensions::Dimensions) = Dimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::Dimensions) = transpose(dimensions)

# this is used to show `dims` for Qobj and QobjEvo
function _get_dims_string(dimensions::Dimensions)
    dims = dimensions_to_dims(dimensions)
    return "($(string(dims[1])), $(string(dims[2])))"
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing

Base.:(==)(dim1::Dimensions, dim2::Dimensions) = (dim1.to == dim2.to) && (dim1.from == dim2.from)

##########################################
# Space

@doc raw"""
    struct Space <: AbstractSpace
        size::Int
    end

A structure that describes a single space with size equals to `size`.
"""
struct Space <: AbstractSpace
    size::Int

    function Space(size::Int)
        (size < 1) && throw(DomainError(size, "The size of `Space` must be positive integer (≥ 1)."))
        return new(size)
    end
end

function Base.show(io::IO, s::Space)
    print(io, "Space($(s.size))")
    return nothing
end

Base.length(s::Space) = 1

get_size(s::Space) = s.size

dimensions_to_dims(s::Space) = SVector{1, Int}(s.size)

##########################################
# TensorSpace

@doc raw"""
    struct TensorSpace{N, T <: NTuple{N, AbstractSpace}} <: AbstractSpace
        spaces::T
    end

A structure that describes a tensor product of `N` spaces, where each space is an element of the tuple `spaces`.
"""
struct TensorSpace{N, T <: NTuple{N, AbstractSpace}} <: AbstractSpace
    spaces::T # a tuple which all elements should be <: AbstractSpace

    function TensorSpace(spaces::T) where {N, T <: NTuple{N, AbstractSpace}}
        if N < 1
            throw(DomainError(N, "The number of spaces in `TensorSpace` must be larger or equal to 1."))
        elseif N == 1
            return spaces[1] # don't wrap with TensorSpace if there is only one space
        else
            return new{N, T}(spaces)
        end
    end
end
TensorSpace(spaces::AbstractSpace) = spaces # don't wrap with TensorSpace if there is only one space
TensorSpace(s::AbstractSpace...) = TensorSpace(s) # this allows convenient function call: TensorSpace(s1, s2, ...)

function Base.show(io::IO, s::TensorSpace)
    print(io, "TensorSpace(")
    join(io, s.spaces, ", ") 
    print(io, ")")
    return nothing
end

Base.length(::TensorSpace{N}) where {N} = N

get_size(s::TensorSpace) = prod(get_size, s.spaces)

Base.kron(s1::TensorSpace, s2::TensorSpace) = TensorSpace(s1.spaces..., s2.spaces...)
Base.kron(s1::TensorSpace, s2::AbstractSpace) = TensorSpace(s1.spaces..., s2)
Base.kron(s1::AbstractSpace, s2::TensorSpace) = TensorSpace(s1, s2.spaces...)
Base.kron(s1::AbstractSpace, s2::AbstractSpace) = TensorSpace(s1, s2)

dimensions_to_dims(s::TensorSpace) = vcat(map(dimensions_to_dims, s.spaces)...)

##########################################
# LiouvilleSpace

@doc raw"""
    LiouvilleSpace{T <: Dimensions} <: AbstractSpace
        op_dims::T
    end

A structure that describes the Liouville space which is equivalent (isometrically isomorphic) to the tensor product of the original `oper` with its dual.
"""
struct LiouvilleSpace{DT <: Dimensions} <: AbstractSpace
    op_dims::DT # original operator dimensions

    function LiouvilleSpace(op_dims::Dimensions{T1, T2}) where {T1 <: AbstractSpace, T2 <: AbstractSpace}
        T1 != T2 && throw(ArgumentError("`LiouvilleSpace` only allows equal `to` and `from` spaces in `Dimensions`."))
        return new{T1}(op_dims)
    end
end

function Base.show(io::IO, s::LiouvilleSpace)
    print(io, "LiouvilleSpace(", s.op_dims, ")")
    return nothing
end

Base.length(::LiouvilleSpace) = 1

get_size(s::LiouvilleSpace) = prod(get_size(s.op_dims)) # get_size(Dimensions.to) × get_size(Dimensions.from)

dimensions_to_dims(s::LiouvilleSpace) = dimensions_to_dims(s.op_dims)
