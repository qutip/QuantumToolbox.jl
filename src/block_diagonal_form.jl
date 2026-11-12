export BlockDiagonalForm, block_diagonal_form

@doc raw"""
    struct BlockDiagonalForm

A type for storing a block-diagonal form of a matrix.

# Fields
- `B`: The block-diagonal matrix. It can be a sparse matrix or a [`QuantumObject`](@ref).
- `P`: The permutation matrix, such that `B = P * A * P'`. It can be a sparse matrix or a [`QuantumObject`](@ref).
- `blocks::AbstractVector`: The blocks of the block-diagonal matrix.
- `block_sizes::AbstractVector`: The sizes of the blocks.
"""
struct BlockDiagonalForm{BT, PT, BLT <: AbstractVector, BST <: AbstractVector}
    B::BT
    P::PT
    blocks::BLT
    block_sizes::BST
end

# indices of each block: connected components of the graph of the nonzero elements of A
_block_indices(A::AbstractSparseMatrix) = connected_components(DiGraph(abs.(A)))
_block_indices(A::AbstractMatrix) = _block_indices(sparse(A))

# permutation matrix that puts the indices of the blocks one after the other
function _block_permutation(A::MT, indices) where {MT <: AbstractMatrix}
    n = LinearAlgebra.checksquare(A)
    return MT(sparse(1:n, reduce(vcat, indices), ones(n), n, n))
end

function block_diagonal_form(A::MT) where {MT <: AbstractSparseMatrix}
    n = LinearAlgebra.checksquare(A)
    idxs_list = _block_indices(A)
    block_sizes = length.(idxs_list)
    P = _block_permutation(A, idxs_list)
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

Return the block-diagonal form of an [`Operator`](@ref) or a [`SuperOperator`](@ref) `A`. This is very useful in the presence of symmetries: see the [Symmetries and block diagonalization](@ref doc:Symmetries) page of the documentation.

The blocks are the connected components of the graph of the nonzero elements of `A`. The ``k``-th block acts on a space ``S_k``, and `A` in the basis of the blocks acts on the direct sum ``S = S_1 \oplus S_2 \oplus \cdots`` (see [`SumSpace`](@ref)).

# Arguments
- `A::QuantumObject`: The quantum object.

# Returns
The [`BlockDiagonalForm`](@ref) `bdf` of `A`, with
- `bdf.B`: `A` in the basis of the blocks, with dims `(S, S)`;
- `bdf.P`: the permutation from the basis of `A` to the basis of the blocks, with dims `(S, dims of A)`, such that `bdf.B = bdf.P * A * bdf.P'`;
- `bdf.blocks`: the blocks, with dims `(Sₖ, Sₖ)`, such that `bdf.B == direct_sum(bdf.blocks...)`;
- `bdf.block_sizes`: the sizes of the blocks.

For an [`Operator`](@ref), each ``S_k`` is a [`Space`](@ref).

For a [`SuperOperator`](@ref) acting on the operators on a Hilbert space ``V`` (e.g., a Liouvillian), the blocks also split ``V`` into sectors ``V_1 \oplus V_2 \oplus \cdots``, and each ``S_k`` is a direct sum of the [`LiouvilleSpace`](@ref)s of some blocks ``\hat{\rho}_{ij}`` of the density matrix (the operators from ``V_j`` to ``V_i``). The vector of a block is made of the vectorized ``\hat{\rho}_{ij}``, one after the other, so that `vec2mat(direct_component(v, i))` returns the ``i``-th of them as a matrix. As direct sums are flattened, ``S`` is the direct sum of all these Liouville spaces, and `direct_sum(S, k => v)` brings a vector `v` of a block back to the space of `bdf.B` if the block starts at the `k`-th summand of ``S``.
"""
function block_diagonal_form(A::QuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}}
    indices = _block_indices(A.data)
    return _block_diagonal_form(A, indices, map(idx -> Space(length(idx)), indices))
end

function block_diagonal_form(L::QuantumObject{SuperOperator, <:Dimensions{<:LiouvilleSpace, <:LiouvilleSpace}})
    isendomorphic(L.dimensions.to.op_dims) || _non_endomorphic_dims_error("operators the SuperOperator acts on", L.dimensions.to.op_dims)
    V = L.dimensions.to.op_dims.to
    n = get_size(V)
    sectors = _block_indices(L.data)

    # 1. C[i, j] is the block of the basis operator |i⟩⟨j| (its index in mat2vec is the linear index of C[i, j])
    C = zeros(Int, n, n)
    for (k, idx) in enumerate(sectors)
        C[idx] .= k
    end

    # 2. Hilbert sectors V₁, V₂, …: the basis states with the same row and the same column of C. C is then constant
    #    on each Vᵢ × Vⱼ, so every block of L is made of whole blocks ρᵢⱼ of the density matrix.
    keys = [(C[i, :], C[:, i]) for i in 1:n]
    hilbert = [findall(==(key), keys) for key in unique(keys)]

    # 3. the blocks (Vᵢ, Vⱼ) of each block of L, with their indices in the order of mat2vec
    pieces = [[(I, J) for I in hilbert for J in hilbert if C[first(I), first(J)] == k] for k in eachindex(sectors)]
    indices = [reduce(vcat, [vec(I .+ n .* (J' .- 1)) for (I, J) in b]) for b in pieces]
    spaces = [SumSpace([LiouvilleSpace(Dimensions(Space(length(I)), Space(length(J)))) for (I, J) in b]) for b in pieces]
    return _block_diagonal_form(L, indices, spaces)
end

# `spaces[k]` is the space of the k-th block, and the space of B is the (flattened) direct sum of them
function _block_diagonal_form(A::QuantumObject, indices, spaces)
    isendomorphic(A.dimensions) || _non_endomorphic_dims_error("QuantumObject", A.dimensions)

    S = SumSpace(reduce(vcat, [collect(_summands(s)) for s in spaces]))
    P = QuantumObject(_block_permutation(A.data, indices), A.type, Dimensions(S, A.dimensions.to))
    B = QuantumObject(P.data * A.data * P.data', A.type, Dimensions(S, S))

    block_sizes = length.(indices)
    offsets = cumsum(vcat(0, block_sizes))
    blocks = map(eachindex(indices)) do k
        r = (offsets[k] + 1):offsets[k + 1]
        return QuantumObject(B.data[r, r], A.type, Dimensions(spaces[k], spaces[k]))
    end
    return BlockDiagonalForm(B, P, blocks, block_sizes)
end
