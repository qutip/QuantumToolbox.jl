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
function direct_sum(S::SumSpace, components::Pair{Int, <:QuantumObject}...)
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
