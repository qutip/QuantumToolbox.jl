using Graphs
import Graphs.Parallel

function bdf(A::SparseMatrixCSC{T,M}) where {T,M}
    n = LinearAlgebra.checksquare(A)

    G = DiGraph(abs.(A) .> 0)
    # G = Graph(A+A')
    v1 = 1
    S = Parallel.bfs_tree(G, v1)
    V = vcat([v1], map(i->dst(i), edges(S)))
    m = length(V)
    m == n && return (spdiagm(ones(n)), A, [n])
    i = 1:n
    V_c = setdiff(i, V)

    block_sizes = [m]
    
    while m < n
        v1 = V_c[1]
        S = Parallel.bfs_tree(G, v1)
        V = vcat(V, [v1], map(i->dst(i), edges(S)))
        V_c = setdiff(i, V)
        m = length(V)
        push!(block_sizes, m-sum(block_sizes))
    end
    
    P = sparse(i, V, ones(n), n, n)
    P, P * A * P', block_sizes
end

function bdf(A::QuantumObject{SparseMatrixCSC{T,M},OpType}) where {T,M,OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}}
    P, A_bd, block_sizes = bdf(A.data)
    P, QuantumObject(A_bd, A.type, A.dims), block_sizes
end

function get_bdf_blocks(A::SparseMatrixCSC{T,M}, block_sizes::Vector{Int}) where {T,M}
    num_blocks = length(block_sizes)
    block_indices = M[1]
    block_list = [A[1:block_sizes[1],1:block_sizes[1]]]
    for i in 2:num_blocks
        idx = sum(view(block_sizes, 1:i-1)) + 1
        push!(block_indices, idx)
        push!(block_list, A[idx:idx-1+block_sizes[i],idx:idx-1+block_sizes[i]])
    end
    block_list, block_indices
end

function get_bdf_blocks(A::QuantumObject{SparseMatrixCSC{T,M},OpType}, block_sizes::Vector{Int}) where {T,M,OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}}
    get_bdf_blocks(A.data, block_sizes)
end