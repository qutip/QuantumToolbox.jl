struct LocalTensorProductOperator{T, M, N, D <: NTuple{N}, I <: NTuple{M}, O <: Tuple{Vararg{AbstractSciMLOperator{T}, M}}, C} <: AbstractSciMLOperator{T}
    dims::D
    indices::I
    ops::O
    cache::C
end

"""
    LocalTensorProductOperator(dims::Tuple, pairs::Pair{Int}...)

Sparse constructor: specify only the non-identity subsystem operators.

# Example
```julia
LocalTensorProductOperator((4, 3, 5), 1 => A, 3 => B)
```
"""
function LocalTensorProductOperator(dims::Tuple, pairs::Pair{Int, <:AbstractSciMLOperator}...)
    if isempty(pairs)
        return LocalTensorProductOperator(dims, (), (), nothing)
    end

    # This keeps them as a Tuple and so we keep the information at compile time
    indices = map(first, pairs)
    ops = map(last, pairs)

    return LocalTensorProductOperator(dims, indices, ops, nothing)
end

Base.size(L::LocalTensorProductOperator) = (prod(L.dims), prod(L.dims))
Base.size(L::LocalTensorProductOperator, n::Int) = size(L)[n]

SciMLOperators.islinear(::LocalTensorProductOperator) = true
SciMLOperators.has_adjoint(L::LocalTensorProductOperator) = all(has_adjoint, L.ops)

function Base.adjoint(L::LocalTensorProductOperator)
    adj_ops = map(adjoint, L.ops)
    return LocalTensorProductOperator(L.dims, L.indices, adj_ops, L.cache)
end

_needs_cache(L::LocalTensorProductOperator{T, M, N}) where {T, M, N} = (M > 1) || (M == 1 && L.indices[1] != N)

SciMLOperators.iscached(L::LocalTensorProductOperator{T, M, N}) where {T, M, N} = !_needs_cache(L) || !isnothing(L.cache)

function cache_operator(L::LocalTensorProductOperator{T, M, N}, v::AbstractVector) where {T, M, N}
    total = prod(L.dims)

    buf = _needs_cache(L) ? (similar(v, total), similar(v, total)) : nothing

    # Also cache internals of each sub-operator
    cached_ops = ntuple(Val(M)) do j
        dk = L.dims[L.indices[j]]
        u = @view(v[1:dk])
        cache_operator(L.ops[j], u)
    end

    return LocalTensorProductOperator(L.dims, L.indices, cached_ops, buf)
end

function LinearAlgebra.mul!(w::AbstractVector, L::LocalTensorProductOperator{T, M, N}, v::AbstractVector) where {T, M, N}
    dims = L.dims
    indices = L.indices
    ops = L.ops

    dims_rev = reverse(dims)

    if M == 0
        copyto!(w, v)
        return w
    end

    # M >= 1
    iscached(L) || throw(ArgumentError("Operator is not cached"))
    buf = L.cache[1]

    current_src = v
    for j in 1:M
        idx = N - indices[j] + 1

        # Ping-pong between buf and w to avoid unnecessary allocations
        current_dst = iseven(j) ? buf : w
        current_buf = iseven(j) ? w : buf

        _apply_single_op_tensor_prod!(current_dst, ops[j], current_src, dims_rev, idx, current_buf)
        current_src = current_dst
    end

    # If M is even, the last write went to buf — copy to w
    if iseven(M)
        copyto!(w, buf)
    end

    return w
end

function _apply_single_op_tensor_prod!(
        dst::AbstractVector, op, src::AbstractVector,
        dims_rev::NTuple{N, Int}, idx::Int, ::Nothing,
    ) where {N}
    dk = dims_rev[idx]
    rest = length(src) ÷ dk

    # If the cache is nothing, we know that the target dimension is already contiguous,
    # so we can skip the permutation step
    mul!(reshape(dst, dk, rest), op, reshape(src, dk, rest))
    return dst
end

function _apply_single_op_tensor_prod!(
        dst::AbstractVector, op, src::AbstractVector,
        dims_rev::NTuple{N, Int}, idx::Int, perm_buf::AbstractVector,
    ) where {N}
    dk = dims_rev[idx]
    rest = length(src) ÷ dk

    if idx == 1
        # Fast path: target dimension is already contiguous
        mul!(reshape(dst, dk, rest), op, reshape(src, dk, rest))
    else
        perm = ntuple(Val(N)) do j
            if j == 1
                return idx
            else
                return (j ≤ idx) ? j - 1 : j
            end
        end
        inv_perm = ntuple(Val(N)) do j
            if j == idx
                return 1
            else
                return (j < idx) ? j + 1 : j
            end
        end
        perm_dims = ntuple(i -> dims_rev[perm[i]], Val(N))

        permutedims!(reshape(dst, perm_dims...), reshape(src, dims_rev...), perm)
        mul!(reshape(perm_buf, dk, rest), op, reshape(dst, dk, rest))
        permutedims!(reshape(dst, dims_rev...), reshape(perm_buf, perm_dims...), inv_perm)
    end
    return dst
end

function LinearAlgebra.mul!(
        w::AbstractVector,
        L::LocalTensorProductOperator{T, M, N},
        v::AbstractVector,
        α,
        β,
    ) where {T, M, N}
    length(w) == prod(L.dims) || throw(DimensionMismatch("output vector has wrong length"))
    length(v) == prod(L.dims) || throw(DimensionMismatch("input vector has wrong length"))

    # Fast exits for scalar coefficients
    if iszero(α)
        if iszero(β)
            fill!(w, zero(eltype(w)))
        elseif !isone(β)
            rmul!(w, β)
        end
        return w
    end

    # Identity operator case: w <- α*v + β*w
    if M == 0
        return axpby!(α, v, β, w)
    end

    # If β == 0, compute w <- L*v first, then scale by α if needed
    if iszero(β)
        mul!(w, L, v)
        if !isone(α)
            rmul!(w, α)
        end
        return w
    end

    # General case:
    # compute tmp = L*v into an internal cached buffer, then w <- α*tmp + β*w
    iscached(L) || throw(ArgumentError("Operator is not cached"))
    buf2 = L.cache[2]

    mul!(buf2, L, v)
    axpby!(α, buf2, β, w)
    return w
end
