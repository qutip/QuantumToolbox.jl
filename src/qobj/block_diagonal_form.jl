export BlockDiagonalForm, block_diagonal_form

@doc raw"""
    struct BlockDiagonalForm

A type for storing a block-diagonal form of a matrix.

# Fields
- `B::DT`: The block-diagonal matrix. It can be a sparse matrix or a [`QuantumObject`](@ref).
- `P::DT`: The permutation matrix. It can be a sparse matrix or a [`QuantumObject`](@ref).
- `blocks::AbstractVector`: The blocks of the block-diagonal matrix.
- `block_sizes::AbstractVector`: The sizes of the blocks.
"""
struct BlockDiagonalForm{DT,BT<:AbstractVector,BST<:AbstractVector}
    B::DT
    P::DT
    blocks::BT
    block_sizes::BST
end

function block_diagonal_form(A::MT) where {MT<:AbstractSparseMatrix}
    n = LinearAlgebra.checksquare(A)

    G = DiGraph(abs.(A))
    idxs_list = connected_components(G)
    block_sizes = length.(idxs_list)

    P = MT(sparse(1:n, reduce(vcat, idxs_list), ones(n), n, n))

    blocks = map(eachindex(idxs_list)) do i
        m = block_sizes[i]
        idxs = idxs_list[i]
        P_i = MT(sparse(1:m, idxs, ones(m), m, n))
        return P_i * A * P_i'
    end

    B = P * A * P'

    return BlockDiagonalForm(B, P, blocks, block_sizes)
end

@doc raw"""
    block_diagonal_form(A::QuantumObject)

Return the block-diagonal form of a [`QuantumObject`](@ref). This is very useful in the presence of symmetries.

# Arguments
- `A::QuantumObject`: The quantum object.

# Returns
The [`BlockDiagonalForm`](@ref) of `A`.
"""
function block_diagonal_form(
    A::QuantumObject{OpType},
) where {OpType<:Union{Operator,SuperOperator}}
    bdf = block_diagonal_form(A.data)
    B = QuantumObject(bdf.B, type = A.type, dims = A.dimensions)
    P = QuantumObject(bdf.P, type = A.type, dims = A.dimensions)
    return BlockDiagonalForm(B, P, bdf.blocks, bdf.block_sizes)
end
