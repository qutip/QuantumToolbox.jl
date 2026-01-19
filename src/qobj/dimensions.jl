#=
This file defines the ProductDimensions structure, which can describe composite Hilbert spaces.
=#

export AbstractDimensions, ProductDimensions
export get_hilbert_size, get_liouville_size

abstract type AbstractDimensions{M, N} end

@doc raw"""
    struct ProductDimensions{M, N, T1<:Tuple, T2<:Tuple} <: AbstractDimensions{M, N}
        to::T1
        from::T2
    end

A structure that embodies the left-hand side (`to`) and right-hand side (`from`) [`AbstractSpace`](@ref) of a quantum object.

For square quantum objects (where `to == from`), such as standard operators, kets, bras, superoperators, etc.,
the `to` and `from` fields will be identical.

For non-square operators, `to` represents the row dimensions and `from` represents the column dimensions.

# Constructors

- `ProductDimensions(to, from)`: Create with explicit `to` and `from` spaces (tuples of `AbstractSpace`)
- `ProductDimensions(dims)`: Create square dimensions where `to == from`
- `ProductDimensions((to_dims, from_dims))`: Create non-square dimensions from a 2-element tuple of integer vectors

# Examples

```jldoctest
julia> ProductDimensions(3)  # Single 3-dimensional Hilbert space
ProductDimensions{1, 1, Tuple{HilbertSpace}, Tuple{HilbertSpace}}((HilbertSpace(3),), (HilbertSpace(3),))

julia> ProductDimensions((2, 3))  # Composite 2⊗3 Hilbert space (square)
ProductDimensions{2, 2, Tuple{HilbertSpace, HilbertSpace}, Tuple{HilbertSpace, HilbertSpace}}((HilbertSpace(2), HilbertSpace(3)), (HilbertSpace(2), HilbertSpace(3)))

julia> ProductDimensions(([2, 3], [4]))  # Non-square: maps from 4-dim to 2⊗3=6-dim
ProductDimensions{2, 1, Tuple{HilbertSpace, HilbertSpace}, Tuple{HilbertSpace}}((HilbertSpace(2), HilbertSpace(3)), (HilbertSpace(4),))
```
"""
struct ProductDimensions{M, N, T1 <: Tuple, T2 <: Tuple} <: AbstractDimensions{M, N}
    to::T1   # space acting on the left (rows)
    from::T2 # space acting on the right (columns)

    # make sure the elements in the tuple are all AbstractSpace
    function ProductDimensions(to::NTuple{M, AbstractSpace}, from::NTuple{N, AbstractSpace}) where {M, N}
        return new{M, N, typeof(to), typeof(from)}(to, from)
    end
end

# Square dimensions constructor from tuple of AbstractSpace
function ProductDimensions(dims::NTuple{N, AbstractSpace}) where {N}
    return ProductDimensions(dims, dims)
end

# Square dimensions from integer tuple/vector
function ProductDimensions(dims::Union{AbstractVector{T}, NTuple{N, T}}) where {T <: Integer, N}
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    spaces = Tuple(HilbertSpace.(dims))
    return ProductDimensions(spaces, spaces)
end

# Square dimensions from single integer
ProductDimensions(dims::Int) = ProductDimensions((HilbertSpace(dims),))

# Square dimensions from single AbstractSpace
ProductDimensions(dims::DimType) where {DimType <: AbstractSpace} = ProductDimensions((dims,))

# Non-square dimensions from 2-element tuple of integer vectors/tuples
function ProductDimensions(dims::Union{AbstractVector{T}, NTuple{2, T}}) where {T <: Union{AbstractVector, NTuple}}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    _non_static_array_warning("dims[1]", dims[1])
    _non_static_array_warning("dims[2]", dims[2])

    L1 = length(dims[1])
    L2 = length(dims[2])
    (L1 > 0) || throw(DomainError(L1, "The length of `dims[1]` must be larger or equal to 1."))
    (L2 > 0) || throw(DomainError(L2, "The length of `dims[2]` must be larger or equal to 1."))

    return ProductDimensions(Tuple(HilbertSpace.(dims[1])), Tuple(HilbertSpace.(dims[2])))
end

# Error for invalid input
ProductDimensions(dims::Any) = throw(
    ArgumentError(
        "The argument dims must be a Tuple or a StaticVector of non-zero length and contain only positive integers.",
    ),
)

# Check if dimensions are square (to == from)
issquare(dimensions::ProductDimensions) = dimensions.to == dimensions.from

_gen_dimensions(dims::AbstractDimensions) = dims
_gen_dimensions(dims::Union{AbstractVector{T}, NTuple{N, T}}) where {T <: Integer, N} = ProductDimensions(dims)
_gen_dimensions(dims::Union{AbstractVector{T}, NTuple{N, T}}) where {T <: Union{AbstractVector, NTuple}, N} =
    ProductDimensions(dims)
_gen_dimensions(dims::Any) = ProductDimensions(dims)

# obtain dims in the type of SVector with integers
dimensions_to_dims(dimensions::NTuple{N, AbstractSpace}) where {N} = vcat(map(dimensions_to_dims, dimensions)...)
function dimensions_to_dims(dimensions::ProductDimensions)
    to_dims = dimensions_to_dims(dimensions.to)
    from_dims = dimensions_to_dims(dimensions.from)
    if dimensions.to == dimensions.from
        return to_dims
    else
        return (to_dims, from_dims)
    end
end

dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

Base.length(::AbstractDimensions{M, N}) where {M, N} = M  # length refers to the number of subsystems in `to`

"""
    get_hilbert_size(dimensions)

Returns the matrix dimensions `(m, n)` of an [`Operator`](@ref) with the given `dimensions`.

Returns `(m, n)` where `m` is the product of the `to` Hilbert space dimensions
and `n` is the product of the `from` Hilbert space dimensions.

If `dimensions` is an `Integer` or a vector/tuple of `Integer`s, it is automatically treated as square [`ProductDimensions`](@ref).
"""
function get_hilbert_size(dimensions::ProductDimensions)
    m = prod(get_hilbert_size, dimensions.to)
    n = prod(get_hilbert_size, dimensions.from)
    return (m, n)
end
get_hilbert_size(dimensions::Union{<:Integer, AbstractVector{<:Integer}, NTuple{N, Integer}}) where {N} =
    get_hilbert_size(ProductDimensions(dimensions))

"""
    get_liouville_size(dimensions)

Returns the matrix dimensions `(m, n)` of a [`SuperOperator`](@ref) with the given `dimensions`.

Returns `(m, n)` where `m` is the product of the `to` Liouville space dimensions
and `n` is the product of the `from` Liouville space dimensions.
(Each Hilbert dimension `d` contributes `d²` to the product.)

If `dimensions` is an `Integer` or a vector/tuple of `Integer`s, it is automatically treated as square [`ProductDimensions`](@ref).
"""
function get_liouville_size(dimensions::ProductDimensions)
    m = prod(get_liouville_size, dimensions.to)
    n = prod(get_liouville_size, dimensions.from)
    return (m, n)
end
get_liouville_size(dimensions::Union{<:Integer, AbstractVector{<:Integer}, NTuple{N, Integer}}) where {N} =
    get_liouville_size(ProductDimensions(dimensions))

Base.transpose(dimensions::ProductDimensions) = ProductDimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::AbstractDimensions) = transpose(dimensions)

# this is used to show `dims` for Qobj and QobjEvo
function _get_dims_string(dimensions::ProductDimensions)
    if issquare(dimensions)
        return string(dimensions_to_dims(dimensions.to))
    else
        dims = dimensions_to_dims(dimensions)
        return "[$(string(dims[1])), $(string(dims[2]))]"
    end
end
_get_dims_string(::Nothing) = "nothing" # for EigsolveResult.dimensions = nothing

Base.:(==)(dim1::ProductDimensions, dim2::ProductDimensions) =
    (dim1.to == dim2.to) && (dim1.from == dim2.from)
