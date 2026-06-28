# ──────────────────────────────────────────────────────────────────────────────
# KroneckerOperator — lazy tensor product of AbstractSciMLOperators
# ──────────────────────────────────────────────────────────────────────────────

"""
    KroneckerOperator{T,D,I,O,C} <: AbstractSciMLOperator{T}

Lazy tensor (Kronecker) product of operators acting on a composite Hilbert space.

Given subsystem dimensions `dims = (d₁, d₂, …, dₙ)`, the operator represents

    O = I_{d₁} ⊗ ⋯ ⊗ Aⱼ ⊗ ⋯ ⊗ I_{dₙ}

where non-identity operators `Aⱼ` are stored only at the subsystem indices
given by `indices`. Identity factors are never materialised.

# Memory layout convention

Physics mode `j` (1 = outermost / leftmost in the Kronecker product,
`n` = innermost / rightmost) maps to Julia array dimension `n - j + 1`
when the flat state vector is reshaped via `reshape(v, reverse(dims)...)`.

*   Physics mode `n` (innermost) → Julia dim 1 (contiguous in memory)
    ⟹  **no permutation needed** (fast path).
*   Physics mode 1 (outermost)  → Julia dim `n` (most strided)
    ⟹  requires `permutedims!`.

# Constructors

    KroneckerOperator(dims, idx₁ => op₁, idx₂ => op₂, …)   # sparse — identities implicit
    KroneckerOperator((op₁, op₂, …, opₙ))                   # full tuple — IdentityOperator marks identity

# Fields

*   `dims::D`    — `NTuple{n,Int}`, dimensions of every subsystem (physics order).
*   `indices::I` — `NTuple{m,Int}`, sorted physics-mode indices with non-identity ops.
*   `ops::O`     — `NTuple{m, <:AbstractSciMLOperator}`, the active operators.
*   `cache::C`   — `Nothing` (default) or `NTuple{3,AbstractArray}` of work buffers
                   allocated by [`cache_operator`](@ref).
"""
struct KroneckerOperator{T, D <: Tuple, I <: Tuple, O <: Tuple, C} <: AbstractSciMLOperator{T}
    dims::D
    indices::I
    ops::O
    cache::C

    function KroneckerOperator(
            dims::D, indices::I, ops::O, cache::C,
        ) where {D <: Tuple, I <: Tuple, O <: Tuple, C}
        # Validate
        n = length(dims)
        m = length(indices)
        length(ops) == m || throw(ArgumentError("length(indices) must equal length(ops)"))
        for j in 1:m
            idx = indices[j]
            (1 <= idx <= n) || throw(ArgumentError("index $idx out of range 1:$n"))
            size(ops[j], 1) == dims[idx] || throw(
                ArgumentError("operator $j has size $(size(ops[j], 1)) but dims[$idx] = $(dims[idx])"),
            )
            size(ops[j], 1) == size(ops[j], 2) || throw(
                ArgumentError("operator $j must be square, got size $(size(ops[j]))"),
            )
        end
        # Check sorted and unique
        for j in 2:m
            indices[j] > indices[j - 1] || throw(
                ArgumentError("indices must be sorted and unique, got $indices"),
            )
        end

        T = isempty(ops) ? Float64 : promote_type(eltype.(ops)...)
        return new{T, D, I, O, C}(dims, indices, ops, cache)
    end
end

# ─── Constructors ────────────────────────────────────────────────────────────

"""
    KroneckerOperator(dims::Tuple, pairs::Pair{Int}...)

Sparse constructor: specify only the non-identity subsystem operators.

# Example
```julia
KroneckerOperator((4, 3, 5), 1 => A, 3 => B)
```
"""
function KroneckerOperator(dims::Tuple, pairs::Pair{Int, <:AbstractSciMLOperator}...)
    if isempty(pairs)
        return KroneckerOperator(dims, (), (), nothing)
    end
    sorted = sort(collect(pairs); by = first)
    indices = Tuple(p.first for p in sorted)
    ops = Tuple(p.second for p in sorted)
    return KroneckerOperator(dims, indices, ops, nothing)
end

"""
    KroneckerOperator(ops::Tuple)

Full-tuple constructor: every subsystem is specified.  `IdentityOperator`
entries are detected and skipped (treated as identity).

# Example
```julia
KroneckerOperator((A, IdentityOperator(3), B))
```
"""
function KroneckerOperator(all_ops::Tuple)
    n = length(all_ops)
    dims = ntuple(j -> size(all_ops[j], 1), n)
    active_idx = Int[]
    active_ops = AbstractSciMLOperator[]
    for j in 1:n
        op = all_ops[j]
        if !(op isa IdentityOperator)
            push!(active_idx, j)
            push!(active_ops, op)
        end
    end
    return KroneckerOperator(dims, Tuple(active_idx), Tuple(active_ops), nothing)
end

# ─── SciMLOperator interface ─────────────────────────────────────────────────

Base.size(L::KroneckerOperator) = (prod(L.dims), prod(L.dims))
Base.size(L::KroneckerOperator, n::Int) = size(L)[n]

islinear(::KroneckerOperator) = true
has_adjoint(L::KroneckerOperator) = all(has_adjoint, L.ops)

function Base.adjoint(L::KroneckerOperator)
    adj_ops = map(adjoint, L.ops)
    return KroneckerOperator(L.dims, L.indices, adj_ops, nothing)
end

# ─── Cache ───────────────────────────────────────────────────────────────────

function cache_operator(L::KroneckerOperator, v::AbstractVecOrMat)
    total = prod(L.dims)
    buf_a = similar(v, total)
    buf_b = similar(v, total)
    buf_c = similar(v, total)
    # Also cache internals of each sub-operator
    cached_ops = _cache_sub_ops(L.ops, L.dims, L.indices, v)
    return KroneckerOperator(L.dims, L.indices, cached_ops, (buf_a, buf_b, buf_c))
end

function _cache_sub_ops(ops::Tuple, dims::Tuple, indices::Tuple, v::AbstractVecOrMat)
    return ntuple(length(ops)) do j
        dk = dims[indices[j]]
        u = similar(v, dk)
        cache_operator(ops[j], u)
    end
end

iscached(L::KroneckerOperator) = L.cache !== nothing

# ─── Core helper: apply single operator to one subsystem ─────────────────────

"""
    _apply_single_op!(dst, op, src, julia_dims, julia_d, perm_buf)

Apply `op` (acting on Julia dimension `julia_d` of the tensor with dimensions
`julia_dims`) to `src`, writing the result to `dst`.

`perm_buf` is scratch space that must not alias `src` or `dst`.
When `julia_d == 1` (innermost physics mode — contiguous in memory),
the permutation is skipped entirely.
"""
function _apply_single_op!(
        dst::AbstractVector, op, src::AbstractVector,
        julia_dims::NTuple{N, Int}, julia_d::Int, perm_buf::AbstractVector,
    ) where {N}
    dk = julia_dims[julia_d]
    rest = length(src) ÷ dk

    if julia_d == 1
        # Fast path: target dimension is already contiguous
        mul!(reshape(dst, dk, rest), op, reshape(src, dk, rest))
    else
        _apply_single_op_permuted!(dst, op, src, julia_dims, Val(julia_d), perm_buf)
    end
    return dst
end

function _apply_single_op_permuted!(
        dst::AbstractVector, op, src::AbstractVector,
        julia_dims::NTuple{N, Int}, ::Val{D}, perm_buf::AbstractVector,
    ) where {N, D}
    dk = julia_dims[D]
    rest = length(src) ÷ dk

    # Build permutation that moves julia_d to position 1
    perm = _make_perm(Val(N), Val(D))
    inv_perm = _make_inv_perm(Val(N), Val(D))
    perm_dims = ntuple(i -> julia_dims[perm[i]], Val(N))

    # 1. Permute src → perm_buf (bring target dim to front)
    permutedims!(reshape(perm_buf, perm_dims...), reshape(src, julia_dims...), perm)

    # 2. Apply operator: perm_buf → dst (both in permuted layout)
    mul!(reshape(dst, dk, rest), op, reshape(perm_buf, dk, rest))

    # 3. Inverse permute: dst (permuted) → perm_buf (original) → dst
    permutedims!(reshape(perm_buf, julia_dims...), reshape(dst, perm_dims...), inv_perm)
    copyto!(dst, perm_buf)

    return dst
end

# ─── Permutation helpers (compile-time via Val) ──────────────────────────────

@generated function _make_perm(::Val{N}, ::Val{D}) where {N, D}
    perm = Vector{Int}(undef, N)
    perm[1] = D
    idx = 2
    for i in 1:N
        if i != D
            perm[idx] = i
            idx += 1
        end
    end
    return :($(Tuple(perm)))
end

@generated function _make_inv_perm(::Val{N}, ::Val{D}) where {N, D}
    # inv of perm = (D, 1, 2, ..., D-1, D+1, ..., N)
    inv_perm = Vector{Int}(undef, N)
    for i in 1:(D - 1)
        inv_perm[i] = i + 1
    end
    inv_perm[D] = 1
    for i in (D + 1):N
        inv_perm[i] = i
    end
    return :($(Tuple(inv_perm)))
end

# ─── 3-arg mul!: w = L * v  (2 buffers) ─────────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVector, L::KroneckerOperator, v::AbstractVector)
    _kron_mul_impl!(w, L, v)
    return w
end

function _kron_mul_impl!(w::AbstractVector, L::KroneckerOperator, v::AbstractVector)
    dims = L.dims
    indices = L.indices
    ops = L.ops
    M = length(indices)
    n = length(dims)
    julia_dims = _reverse_tuple(dims)

    if M == 0
        copyto!(w, v)
        return w
    end

    # Get or allocate buffers
    if L.cache !== nothing
        buf_a, buf_b, _ = L.cache
    else
        buf_a = similar(v)
        buf_b = similar(v)
    end

    current_src = v
    for j in 1:M
        julia_d = n - indices[j] + 1   # physics mode → Julia dimension
        is_last = (j == M)

        if M == 1
            # Single operator: v → w directly, use buf_a as perm scratch
            _apply_single_op!(w, ops[j], v, julia_dims, julia_d, buf_a)
            return w
        end

        if is_last
            current_dst = w
            # Use whichever buffer is NOT current_src as perm scratch
            perm_buf = (current_src === buf_a) ? buf_b : buf_a
        else
            # Ping-pong between buf_a and buf_b
            if j == 1
                current_dst = buf_a
            else
                current_dst = (current_src === buf_a) ? buf_b : buf_a
            end
            # w is not yet written → safe to use as perm scratch
            perm_buf = w
        end

        _apply_single_op!(current_dst, ops[j], current_src, julia_dims, julia_d, perm_buf)
        current_src = current_dst
    end

    return w
end

# ─── 5-arg mul!: w = α * L * v + β * w ──────────────────────────────────────

function LinearAlgebra.mul!(
        w::AbstractVector, L::KroneckerOperator, v::AbstractVector, α, β,
    )
    M = length(L.indices)

    if M == 0
        # All identity: w = α * v + β * w
        if iszero(β)
            copyto!(w, v)
            lmul!(α, w)
        else
            axpby!(α, v, β, w)
        end
        return w
    end

    if iszero(β)
        # w = α * L * v  →  compute L*v into w, then scale
        _kron_mul_impl!(w, L, v)
        lmul!(α, w)
    else
        # Need a 3rd buffer to hold L*v without clobbering w
        if L.cache !== nothing
            five_arg_buf = L.cache[3]
        else
            five_arg_buf = similar(v)
        end
        _kron_mul_impl!(five_arg_buf, L, v)
        axpby!(α, five_arg_buf, β, w)
    end

    return w
end

# ─── Tuple reversal helper ──────────────────────────────────────────────────

@generated function _reverse_tuple(t::NTuple{N, Int}) where {N}
    ex = Expr(:tuple, [:(t[$i]) for i in N:-1:1]...)
    return ex
end

# ─── concretize: materialise the full Kronecker product ──────────────────────

function concretize(L::KroneckerOperator{T}) where {T}
    dims = L.dims
    n = length(dims)
    indices_set = Set(L.indices)

    # Build list of matrices (one per subsystem)
    mats = Vector{AbstractMatrix}(undef, n)
    op_idx = 1
    for j in 1:n
        if j in indices_set
            mats[j] = concretize(L.ops[op_idx])
            op_idx += 1
        else
            d = dims[j]
            mats[j] = sparse(one(T) * I, d, d)
        end
    end

    # Compute kron(mats[1], mats[2], …, mats[n])
    result = mats[1]
    for j in 2:n
        result = kron(result, mats[j])
    end
    return result
end

# ─── Base.kron overloads for KroneckerOperator ───────────────────────────────
# NOTE: We do NOT overload kron(::AbstractSciMLOperator, ::AbstractSciMLOperator)
# because SciMLOperators already defines that → TensorProductOperator.
# Users create KroneckerOperator explicitly via constructors.

# Nested KroneckerOperator flattening: kron(K, op) and kron(op, K)
function Base.kron(A::KroneckerOperator, B::AbstractSciMLOperator)
    n_a = length(A.dims)
    new_dims = (A.dims..., size(B, 1))
    if B isa IdentityOperator
        new_indices = A.indices
        new_ops = A.ops
    else
        new_indices = (A.indices..., n_a + 1)
        new_ops = (A.ops..., B)
    end
    return KroneckerOperator(new_dims, new_indices, new_ops, nothing)
end

function Base.kron(A::AbstractSciMLOperator, B::KroneckerOperator)
    new_dims = (size(A, 1), B.dims...)
    shifted_indices = map(i -> i + 1, B.indices)
    if A isa IdentityOperator
        new_indices = shifted_indices
        new_ops = B.ops
    else
        new_indices = (1, shifted_indices...)
        new_ops = (A, B.ops...)
    end
    return KroneckerOperator(new_dims, new_indices, new_ops, nothing)
end

function Base.kron(A::KroneckerOperator, B::KroneckerOperator)
    n_a = length(A.dims)
    new_dims = (A.dims..., B.dims...)
    shifted_b_indices = map(i -> i + n_a, B.indices)
    new_indices = (A.indices..., shifted_b_indices...)
    new_ops = (A.ops..., B.ops...)
    return KroneckerOperator(new_dims, new_indices, new_ops, nothing)
end

# ─── show ────────────────────────────────────────────────────────────────────

function Base.show(io::IO, L::KroneckerOperator{T}) where {T}
    n = length(L.dims)
    m = length(L.indices)
    return print(io, "KroneckerOperator{$T}(dims=$(L.dims), $m active of $n subsystems)")
end
