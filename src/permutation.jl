export bdf, get_bdf_blocks

function bdf(A::SparseMatrixCSC{T,M}) where {T,M}
    n = LinearAlgebra.checksquare(A)

    G = DiGraph(abs.(A) .> 0)
    idxs = connected_components(G)
    P = sparse(1:n, reduce(vcat, idxs), ones(n), n, n)
    block_sizes = map(length, idxs)

    return P, P * A * P', block_sizes
end

function bdf(
    A::QuantumObject{SparseMatrixCSC{T,M},OpType},
) where {T,M,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    P, A_bd, block_sizes = bdf(A.data)
    return P, QuantumObject(A_bd, A.type, A.dims), block_sizes
end

function get_bdf_blocks(A::SparseMatrixCSC{T,M}, block_sizes::Vector{Int}) where {T,M}
    num_blocks = length(block_sizes)
    block_indices = M[1]
    block_list = [A[1:block_sizes[1], 1:block_sizes[1]]]
    for i in 2:num_blocks
        idx = sum(view(block_sizes, 1:i-1)) + 1
        push!(block_indices, idx)
        push!(block_list, A[idx:idx-1+block_sizes[i], idx:idx-1+block_sizes[i]])
    end
    return block_list, block_indices
end

function get_bdf_blocks(
    A::QuantumObject{SparseMatrixCSC{T,M},OpType},
    block_sizes::Vector{Int},
) where {T,M,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    return get_bdf_blocks(A.data, block_sizes)
end
