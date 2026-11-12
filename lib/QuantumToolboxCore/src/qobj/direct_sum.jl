#=
Direct sums of quantum objects, which live in (or act on) a `SumSpace`.
=#

export direct_sum, direct_component

@doc raw"""
    direct_sum(A₁, A₂, ...)
    A₁ ⊕ A₂

Return the direct sum of spaces, or of [`QuantumObject`](@ref)s of the same type:
- spaces: the [`SumSpace`](@ref) of them;
- [`Ket`](@ref)s (or [`OperatorKet`](@ref)s): the vector made of the given ones, one after the other;
- [`Operator`](@ref)s (or [`SuperOperator`](@ref)s): the block-diagonal matrix made of the given ones.

The direct sum is associative, so a direct sum of direct sums is flattened, as [`kron`](@ref) does for the tensor products: `A₁ ⊕ A₂ ⊕ A₃` is the same as `direct_sum(A₁, A₂, A₃)`.

`⊕` (`\oplus` + `<tab>`) is a synonym of `direct_sum`. See also [`direct_component`](@ref).

# Examples

```jldoctest
julia> H = sigmaz() ⊕ qeye(3) ⊕ sigmax();

julia> H.dimensions
Dimensions(SumSpace(Space(2), Space(3), Space(2)), SumSpace(Space(2), Space(3), Space(2)))

julia> ψ = basis(2, 0) ⊕ basis(3, 1);

julia> direct_component(ψ, 2) == basis(3, 1)
true
```
"""
direct_sum(spaces::AbstractSpace...) = SumSpace(mapreduce(s -> Tuple(_summands(s)), (s1, s2) -> (s1..., s2...), spaces))

function direct_sum(As::QuantumObject{OpType}...) where {OpType <: Union{Operator, SuperOperator}}
    to = direct_sum(map(A -> A.dimensions.to, As)...)
    from = direct_sum(map(A -> A.dimensions.from, As)...)
    return QuantumObject(blockdiag(map(A -> sparse(A.data), As)...), first(As).type, Dimensions(to, from))
end

function direct_sum(vs::QuantumObject{VType}...) where {VType <: Union{Ket, OperatorKet}}
    S = direct_sum(map(v -> v.dimensions.to, vs)...)
    return QuantumObject(vcat(map(v -> v.data, vs)...), first(vs).type, Dimensions(S, Space(1)))
end

# the summands of a space (a space that is not a direct sum is its only summand)
_summands(s::SumSpace) = s.spaces
_summands(s::AbstractSpace) = (s,)

@doc raw"""
    direct_sum(S::SumSpace, k => v, ...)

Return the vector of the [`SumSpace`](@ref) `S` that is equal to the [`Ket`](@ref) (or [`OperatorKet`](@ref)) `v` from the `k`-th component on, and which is zero in the components that are not given.

The space of `v` can be the `k`-th summand of `S`, or the direct sum of some consecutive summands of `S` from the `k`-th on. For example, `v` can be an eigenvector of a block returned by [`block_diagonal_form`](@ref QuantumToolbox.block_diagonal_form).
"""
# (at least one component, otherwise `direct_sum(S)` would be the direct sum of a single space)
function direct_sum(S::SumSpace, c1::Pair{Int, <:QuantumObject}, cs::Pair{Int, <:QuantumObject}...)
    components = (c1, cs...)
    data = zeros(promote_type(map(c -> eltype(c.second), components)...), get_size(S))
    for (k, v) in components
        _check_component(S, k, v.dimensions.to)
        offset = first(_summand_range(S, k)) - 1
        data[(offset + 1):(offset + length(v.data))] .= v.data
    end
    return QuantumObject(data, last(first(components)).type, Dimensions(S, Space(1)))
end

function _check_component(S::SumSpace, k::Int, s::AbstractSpace)
    (S.spaces[k] == s) && return nothing

    # a direct sum of consecutive summands of S
    run = _summands(s)
    r = k:(k + length(run) - 1)
    ((last(r) <= length(S.spaces)) && all(i -> S.spaces[r[i]] == run[i], eachindex(run))) ||
        throw(DimensionMismatch("The dims of the component $k do not match the summands of the SumSpace."))
    return nothing
end

@doc raw"""
    direct_sum(blocks::AbstractMatrix)

Return the [`Operator`](@ref) (or [`SuperOperator`](@ref)) made of a matrix of blocks, where `nothing` stands for a zero block, e.g. `direct_sum([A B; nothing D])`.

The blocks in the same row must have the same `to` space, and the blocks in the same column the same `from` space, so each row and each column needs at least one block (otherwise, see the next method). The result acts from the direct sum of the `from` spaces of the columns to the direct sum of the `to` spaces of the rows, and `direct_component(result, i, j)` is the block `(i, j)`.

# Examples

```jldoctest
julia> a = destroy(4);

julia> direct_sum([qeye(2) nothing; nothing a]) == qeye(2) ⊕ a
true

julia> B = direct_sum([sigmax() basis(2, 0) * basis(4, 1)'; nothing a]);

julia> direct_component(B, 1, 2) == basis(2, 0) * basis(4, 1)'
true
```
"""
function direct_sum(blocks::AbstractMatrix)
    all(b -> isnothing(b) || b isa QuantumObject{<:Union{Operator, SuperOperator}}, blocks) ||
        throw(ArgumentError("The blocks must be Operators, SuperOperators, or `nothing`."))
    to = [_line_space(blocks[i, :], :to, "row $i") for i in axes(blocks, 1)]
    from = [_line_space(blocks[:, j], :from, "column $j") for j in axes(blocks, 2)]
    S_to, S_from = direct_sum(to...), direct_sum(from...)
    return direct_sum(S_to, S_from, (Tuple(idx) => b for (idx, b) in pairs(IndexCartesian(), blocks) if !isnothing(b))...)
end

# the `to` (or `from`) space of the blocks in a row (or column), which must be the same
function _line_space(line, side::Symbol, name::String)
    spaces = [getfield(b.dimensions, side) for b in line if !isnothing(b)]
    isempty(spaces) && throw(ArgumentError("The $name has only zero blocks, so its dims are unknown."))
    allequal(spaces) || throw(DimensionMismatch("The blocks in the $name have different dims."))
    return first(spaces)
end

@doc raw"""
    direct_sum(S_to::SumSpace, S_from::SumSpace, (i, j) => A, ...)

Return the [`Operator`](@ref) (or [`SuperOperator`](@ref)) from `S_from` to `S_to` whose block `(i, j)` is `A`, and which is zero in the blocks that are not given. As for the vectors, the dims of `A` can also span several consecutive summands of `S_to` and `S_from`, starting from the `i`-th and the `j`-th.
"""
# (at least one block, otherwise `direct_sum(S_to, S_from)` would be the direct sum of two spaces)
function direct_sum(S_to::SumSpace, S_from::SumSpace, b1::Pair{<:NTuple{2, Int}, <:QuantumObject}, bs::Pair{<:NTuple{2, Int}, <:QuantumObject}...)
    blocks = (b1, bs...)
    I, J, V = Int[], Int[], promote_type(map(b -> eltype(b.second), blocks)...)[]
    for ((i, j), A) in blocks
        _check_component(S_to, i, A.dimensions.to)
        _check_component(S_from, j, A.dimensions.from)
        Ai, Aj, Av = findnz(sparse(A.data))
        append!(I, Ai .+ (first(_summand_range(S_to, i)) - 1))
        append!(J, Aj .+ (first(_summand_range(S_from, j)) - 1))
        append!(V, Av)
    end
    type = last(first(blocks)).type
    all(b -> b.second.type == type, blocks) || throw(ArgumentError("The blocks must have the same type."))
    return QuantumObject(sparse(I, J, V, get_size(S_to), get_size(S_from)), type, Dimensions(S_to, S_from))
end

@doc raw"""
    direct_component(v::QuantumObject, k)
    direct_component(A::QuantumObject, i, j)

Return the `k`-th component of a [`Ket`](@ref) (or [`OperatorKet`](@ref)) `v` of a [`SumSpace`](@ref), or the block `(i, j)` of an [`Operator`](@ref) (or [`SuperOperator`](@ref)) `A` between two [`SumSpace`](@ref)s. The result has the dims of the corresponding summands.

The indices can also be ranges of consecutive summands, e.g. `direct_component(A, 1:2, 1:2)`, in which case the result acts on the direct sum of these summands.

See also [`direct_sum`](@ref).
"""
function direct_component(
        v::QuantumObject{VType, <:Dimensions{<:SumSpace}},
        k::Union{Int, UnitRange{Int}},
    ) where {VType <: Union{Ket, OperatorKet}}
    S = v.dimensions.to
    return QuantumObject(v.data[_summand_range(S, k)], v.type, Dimensions(_summand_space(S, k), Space(1)))
end

function direct_component(
        A::QuantumObject{OpType, <:Dimensions{<:SumSpace, <:SumSpace}},
        i::Union{Int, UnitRange{Int}},
        j::Union{Int, UnitRange{Int}},
    ) where {OpType <: Union{Operator, SuperOperator}}
    to, from = A.dimensions.to, A.dimensions.from
    data = A.data[_summand_range(to, i), _summand_range(from, j)]
    return QuantumObject(data, A.type, Dimensions(_summand_space(to, i), _summand_space(from, j)))
end

# the space of the k-th summand, or the direct sum of the summands in the range k
_summand_space(S::SumSpace, k::Int) = S.spaces[k]
_summand_space(S::SumSpace, k::UnitRange{Int}) = SumSpace(S.spaces[k])
