#=
This file defines the Dimensions structure and also the following space structures:
    - Space
    - TensorSpace
    - LiouvilleSpace
    - SumSpace
=#

export Dimensions, get_size
export AbstractSpace, Space, TensorSpace, SumSpace
export AbstractSuperSpace, LiouvilleSpace

@doc raw"""
    abstract type AbstractSpace

Abstract type for all space structures.
"""
abstract type AbstractSpace end

@doc raw"""
    abstract type AbstractSuperSpace

Abstract type of space structures for different SuperOperator representation.
"""
abstract type AbstractSuperSpace <: AbstractSpace end

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

- `Dimensions(to, from)`: Create `Dimensions` with explicit `to` and `from` [`AbstractSpace`](@ref)
- `Dimensions(dims)`: Create `Dimensions` where `to == from`
- `Dimensions((to_dims, from_dims))`: Create `Dimensions` from a 2-element tuple of integer tuple/vectors

# Printing

The `dims` of a quantum object are printed as `(to, from)`, where:
- a [`Space`](@ref) is printed as its size, e.g. `2`;
- a [`TensorSpace`](@ref) is printed with `⊗`, e.g. `2⊗3`;
- a [`LiouvilleSpace`](@ref) of the operators with dims `(to, from)`, i.e., of the operators that take a ket of `from` to a ket of `to`, is printed as `ℒ(to ← from)`, e.g. `ℒ(2⊗3 ← 2⊗3)`;
- a [`SumSpace`](@ref) is printed with `⊕`, e.g. `2 ⊕ 3`.

For example, a [`Ket`](@ref) of two qubits has dims `(2⊗2, 1)`, an [`Operator`](@ref) acting on them `(2⊗2, 2⊗2)`, and a [`SuperOperator`](@ref) acting on these operators `(ℒ(2⊗2 ← 2⊗2), ℒ(2⊗2 ← 2⊗2))`. The integer form of the dims is given by `A.dims`.

# Examples

```jldoctest
julia> Dimensions(3)  # Single 3-dimensional Hilbert space
Dimensions(Space(3), Space(3))

julia> Dimensions((2, 3))  # Composite 2⊗3 Hilbert space (square)
Dimensions(TensorSpace(Space(2), Space(3)), TensorSpace(Space(2), Space(3)))

julia> Dimensions(((2, 3), (4,)))  # Non-square: maps from 4-dim to 2⊗3=6-dim
Dimensions(TensorSpace(Space(2), Space(3)), Space(4))
```
"""
struct Dimensions{T1 <: AbstractSpace, T2 <: AbstractSpace}
    to::T1   # space acting on the left (rows)
    from::T2 # space acting on the right (columns)
end

# by defining alias type, it is easier to represent the nested dims structure:
const DimsListType{T1, T2} = Tuple{<:AbstractVecOrTuple{T1}, <:AbstractVecOrTuple{T2}}

function _list_to_tensor_space(dims::AbstractVecOrTuple{T}, argname::String) where {T <: Integer}
    _non_static_array_warning(argname, dims)
    return TensorSpace(Tuple(Space.(dims)))
end

# endomorphic dimensions from integer tuple/vector
function Dimensions(dims::AbstractVecOrTuple{T}) where {T <: Integer}
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))

    return Dimensions(_list_to_tensor_space(dims, "dims"))
end

# endomorphic dimensions from single integer
Dimensions(dims::Int) = Dimensions(Space(dims))

# endomorphic dimensions from single AbstractSpace
Dimensions(dims::AbstractSpace) = Dimensions(dims, dims)

# Non-endomorphic dimensions from 2-element tuple of integer vectors/tuples
function Dimensions(dims::DimsListType{T, T}) where {T <: Integer}
    (length(dims) != 2) && throw(ArgumentError("Invalid dims = $dims"))

    L1 = length(dims[1])
    L2 = length(dims[2])
    (L1 > 0) || throw(DomainError(L1, "The length of `dims[1]` must be larger or equal to 1."))
    (L2 > 0) || throw(DomainError(L2, "The length of `dims[2]` must be larger or equal to 1."))

    return Dimensions(_list_to_tensor_space(dims[1], "dims[1]"), _list_to_tensor_space(dims[2], "dims[2]"))
end

# LiouvilleSpace dimensions for OperatorKet/OperatorBra/SuperOperator from 3-level nested tuple/vector of integers
Dimensions(dims::DimsListType{T1, T2}) where {T1 <: AbstractVecOrTuple, T2 <: Integer} = Dimensions(LiouvilleSpace(Dimensions(dims[1])), _list_to_tensor_space(dims[2], "dims[2]"))
Dimensions(dims::DimsListType{T1, T2}) where {T1 <: Integer, T2 <: AbstractVecOrTuple} = Dimensions(_list_to_tensor_space(dims[1], "dims[1]"), LiouvilleSpace(Dimensions(dims[2])))
Dimensions(dims::DimsListType{T1, T2}) where {T1 <: AbstractVecOrTuple, T2 <: AbstractVecOrTuple} = Dimensions(LiouvilleSpace(Dimensions(dims[1])), LiouvilleSpace(Dimensions(dims[2])))

# Error for invalid input
Dimensions(dims::Any) = throw(
    ArgumentError(
        """
        \nThe argument `dims` must be one of the following:
          • a positive integer (single-system Hilbert space),
          • a `Tuple` or `AbstractVector` of positive integers (multipartite Hilbert space),
          • a `AbstractSpace`,
          • a 2-element `Tuple` combining the above, or
          • directly specify the `Dimensions` structure.
        """
    ),
)

function Base.show(io::IO, d::Dimensions)
    print(io, "Dimensions($(d.to), $(d.from))")
    return nothing
end

# Check if Dimensions is endomorphic (to == from)
isendomorphic(dimensions::Dimensions) = dimensions.to == dimensions.from
_non_endomorphic_dims_error(op_description::String, dimensions::Dimensions) = throw(ArgumentError("The $op_description must be endomorphic (same `to` and `from`), got: dims = $(_get_dims_string(dimensions))"))

# obtain dims in the type of SVector with integers
function dimensions_to_dims(dimensions::Dimensions)
    to_dims = dimensions_to_dims(dimensions.to)
    from_dims = dimensions_to_dims(dimensions.from)
    return (to_dims, from_dims)
end

dimensions_to_dims(::Nothing) = nothing # for EigsolveResult.dimensions = nothing

@doc raw"""
    get_size(dimensions::Dimensions)

Returns the matrix dimensions `(m, n)` of a given [`Dimensions`](@ref).

Returns `(m, n)` where `m` is the product of the `dimensions.to`, and `n` is the product of the `dimensions.from`.

If `dimensions` is an `Integer` or a vector/tuple of `Integer`s, it is automatically treated as `Dimensions(dimensions, dimensions)`.
"""
get_size(dimensions::Dimensions) = (get_size(dimensions.to), get_size(dimensions.from))
get_size(dimensions::Union{T, AbstractVecOrTuple{T}}) where {T <: Integer} = get_size(Dimensions(dimensions))

Base.transpose(dimensions::Dimensions) = Dimensions(dimensions.from, dimensions.to) # switch `to` and `from`
Base.adjoint(dimensions::Dimensions) = transpose(dimensions)

# this is used to show `dims` for Qobj and QobjEvo, in the readable form of `_space_string`
_get_dims_string(dimensions::Dimensions) = "($(_space_string(dimensions.to)), $(_space_string(dimensions.from)))"
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

# allow iteration over subspaces in TensorSpace
Base.iterate(s::TensorSpace) = iterate(s.spaces)
Base.iterate(s::TensorSpace, state) = iterate(s.spaces, state)

Base.getindex(s::TensorSpace, idx::Union{T, AbstractVector{T}}) where {T <: Integer} = TensorSpace(s.spaces[idx]) # indexing for `permute`

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
    LiouvilleSpace{T <: Dimensions} <: AbstractSuperSpace
        op_dims::T
    end

A structure that describes the Liouville space which is equivalent (isometrically isomorphic) to the tensor product of the original `oper` with its dual.
"""
struct LiouvilleSpace{DT <: Dimensions} <: AbstractSuperSpace
    op_dims::DT # original operator dimensions
end

function Base.show(io::IO, s::LiouvilleSpace)
    print(io, "LiouvilleSpace(", s.op_dims, ")")
    return nothing
end

Base.length(::LiouvilleSpace) = 1

# we somehow need to define :(==) for LiouvilleSpace, otherwise, it doesn't work with EnrSpace
Base.:(==)(s1::LiouvilleSpace, s2::LiouvilleSpace) = s1.op_dims == s2.op_dims

get_size(s::LiouvilleSpace) = prod(get_size(s.op_dims)) # get_size(Dimensions.to) × get_size(Dimensions.from)

dimensions_to_dims(s::LiouvilleSpace) = dimensions_to_dims(s.op_dims)

##########################################
# SumSpace

@doc raw"""
    struct SumSpace{T} <: AbstractSpace
        spaces::T
    end

A structure that describes the direct sum ``\mathcal{S}_1 \oplus \mathcal{S}_2 \oplus \cdots`` of the spaces in `spaces` (a `Tuple` or an `AbstractVector`).

Its dimension is the sum of the dimensions of the summands, a vector of it is the concatenation of one vector of each summand, and a matrix acting on it is a block matrix with one block for each pair of summands. Unlike [`TensorSpace`](@ref), a single summand stays a `SumSpace`, so that a block-diagonal matrix with a single block keeps the same structure.

A `SumSpace` is also created by [`direct_sum`](@ref) (or `⊕`) of spaces, which flattens the direct sums of direct sums. It is used by [`direct_sum`](@ref) of quantum objects and by [`block_diagonal_form`](@ref QuantumToolbox.block_diagonal_form).

In the `dims` of a quantum object, a `SumSpace` is printed with `⊕`. For example, `ℒ(2⊗2 ← 2⊗2) ⊕ ℒ(3 ← 3)` is the direct sum of the Liouville spaces of the operators on two qubits and on a three-level system (see [`Dimensions`](@ref)).

# Examples

```jldoctest
julia> S = Space(2) ⊕ Space(3)
SumSpace(Space(2), Space(3))

julia> get_size(S)
5
```
"""
struct SumSpace{T <: Union{Tuple{Vararg{AbstractSpace}}, AbstractVector{<:AbstractSpace}}} <: AbstractSpace
    spaces::T # the summands
end
SumSpace(spaces::AbstractSpace...) = SumSpace(spaces)

function Base.show(io::IO, s::SumSpace)
    print(io, "SumSpace(")
    join(io, s.spaces, ", ")
    print(io, ")")
    return nothing
end

# compare the summands, whether they are stored in a Tuple or in a Vector
Base.:(==)(s1::SumSpace, s2::SumSpace) = (length(s1.spaces) == length(s2.spaces)) && all(splat(==), zip(s1.spaces, s2.spaces))

get_size(s::SumSpace) = sum(get_size, s.spaces)

# the integer form of each summand
dimensions_to_dims(s::SumSpace) = map(dimensions_to_dims, s.spaces)

# the indices of the k-th summand
function _summand_range(s::SumSpace, k::Int)
    offset = sum(get_size, s.spaces[1:(k - 1)]; init = 0)
    return (offset + 1):(offset + get_size(s.spaces[k]))
end

##########################################
# Readable form of the spaces, used to show the `dims` of quantum objects

_space_string(s::AbstractSpace) = sprint(show, s) # e.g. EnrSpace
_space_string(s::Space) = string(s.size)
_space_string(s::TensorSpace) = join(map(_factor_string, s.spaces), "⊗")
_space_string(s::LiouvilleSpace) = "ℒ($(_space_string(s.op_dims.to)) ← $(_space_string(s.op_dims.from)))"
_space_string(s::SumSpace) = join(map(_factor_string, s.spaces), " ⊕ ")

# a space inside a tensor product or a direct sum: direct sums get parentheses
_factor_string(s::AbstractSpace) = _space_string(s)
_factor_string(s::SumSpace) = (length(s.spaces) == 1) ? _space_string(s) : "($(_space_string(s)))"
