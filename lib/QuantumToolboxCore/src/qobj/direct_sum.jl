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

!!! note
    `A₁ ⊕ A₂ ⊕ A₃` is the nested direct sum `(A₁ ⊕ A₂) ⊕ A₃`. Use `direct_sum(A₁, A₂, A₃)` for a direct sum with three summands.

`⊕` (`\oplus` + `<tab>`) is a synonym of `direct_sum`. See also [`direct_component`](@ref).

# Examples

```jldoctest
julia> H = direct_sum(sigmaz(), qeye(3));

julia> H.dimensions
Dimensions(SumSpace(Space(2), Space(3)), SumSpace(Space(2), Space(3)))

julia> ψ = basis(2, 0) ⊕ basis(3, 1);

julia> direct_component(ψ, 2) == basis(3, 1)
true
```
"""
direct_sum(s1::AbstractSpace, s2::AbstractSpace...) = SumSpace(s1, s2...)

function direct_sum(As::QuantumObject{OpType}...) where {OpType <: Union{Operator, SuperOperator}}
    to = SumSpace(map(A -> A.dimensions.to, As))
    from = SumSpace(map(A -> A.dimensions.from, As))
    return QuantumObject(blockdiag(map(A -> sparse(A.data), As)...), first(As).type, Dimensions(to, from))
end

direct_sum(vs::QuantumObject{VType}...) where {VType <: Union{Ket, OperatorKet}} =
    direct_sum(SumSpace(map(v -> v.dimensions.to, vs)), pairs(vs)...)

@doc raw"""
    direct_sum(S::SumSpace, k => v, ...)

Return the vector of the [`SumSpace`](@ref) `S` whose `k`-th component is the [`Ket`](@ref) (or [`OperatorKet`](@ref)) `v`, and which is zero in the components that are not given.
"""
function direct_sum(S::SumSpace, components::Pair{Int, <:QuantumObject}...)
    data = zeros(promote_type(map(c -> eltype(c.second), components)...), get_size(S))
    for (k, v) in components
        (S.spaces[k] == v.dimensions.to) || throw(DimensionMismatch("The dims of the component $k do not match the summand $k of the SumSpace."))
        data[_summand_range(S, k)] .= v.data
    end
    return QuantumObject(data, last(first(components)).type, Dimensions(S, Space(1)))
end

@doc raw"""
    direct_component(v::QuantumObject, k::Int)
    direct_component(A::QuantumObject, i::Int, j::Int)

Return the `k`-th component of a [`Ket`](@ref) (or [`OperatorKet`](@ref)) `v` of a [`SumSpace`](@ref), or the block `(i, j)` of an [`Operator`](@ref) (or [`SuperOperator`](@ref)) `A` between two [`SumSpace`](@ref)s. The result has the dims of the corresponding summands.

See also [`direct_sum`](@ref).
"""
function direct_component(v::QuantumObject{VType, <:Dimensions{<:SumSpace}}, k::Int) where {VType <: Union{Ket, OperatorKet}}
    S = v.dimensions.to
    return QuantumObject(v.data[_summand_range(S, k)], v.type, Dimensions(S.spaces[k], Space(1)))
end

function direct_component(
        A::QuantumObject{OpType, <:Dimensions{<:SumSpace, <:SumSpace}},
        i::Int,
        j::Int,
    ) where {OpType <: Union{Operator, SuperOperator}}
    to, from = A.dimensions.to, A.dimensions.from
    data = A.data[_summand_range(to, i), _summand_range(from, j)]
    return QuantumObject(data, A.type, Dimensions(to.spaces[i], from.spaces[j]))
end
