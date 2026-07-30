using QuantumToolboxOperators
using LinearAlgebra
using SparseArrays
using SciMLOperators
using SciMLOperators: AdjointOperator
using Test

println("="^60)
println("QuantumToolboxOperators — Informal Tests")
println("="^60)

# ── Helper: compare matrix-free mul! against concretize'd sparse ─────────────
function test_mul(L, label; N_trials = 5)
    A = concretize(L)
    N = size(L, 1)
    ok = true
    for _ in 1:N_trials
        v = randn(ComplexF64, N)
        w1 = similar(v)
        mul!(w1, L, v)
        w2 = A * v
        if !isapprox(w1, w2; atol = 1.0e-10)
            ok = false
            @error "$label: 3-arg mul! mismatch" norm(w1 - w2)
        end
    end
    ok && println("  ✓ $label — 3-arg mul! matches concretize")
    return ok
end

function test_mul5(L, label; N_trials = 5)
    A = concretize(L)
    N = size(L, 1)
    ok = true
    for _ in 1:N_trials
        v = randn(ComplexF64, N)
        w = randn(ComplexF64, N)
        w_orig = copy(w)
        α, β = randn(ComplexF64), randn(ComplexF64)
        mul!(w, L, v, α, β)
        w_expected = α * (A * v) + β * w_orig
        if !isapprox(w, w_expected; atol = 1.0e-10)
            ok = false
            @error "$label: 5-arg mul! mismatch" norm(w - w_expected)
        end
    end
    ok && println("  ✓ $label — 5-arg mul! matches concretize")
    return ok
end

# ── Sizes to test ────────────────────────────────────────────────────────────
for N in [2, 5, 10, 20, 50]
    println("\n--- N = $N ---")

    a_op = DestroyOperator(N)
    ad_op = adjoint(a_op)
    n_op0 = NumberOperator(N)
    n_op1 = NumberOperator(N; shift = 1)

    # DestroyOperator
    test_mul(a_op, "DestroyOperator")
    test_mul5(a_op, "DestroyOperator")

    # AdjointOperator (creation)
    test_mul(ad_op, "adjoint(DestroyOperator)")
    test_mul5(ad_op, "adjoint(DestroyOperator)")

    # NumberOperator (shift=0)
    test_mul(n_op0, "NumberOperator(shift=0)")
    test_mul5(n_op0, "NumberOperator(shift=0)")

    # NumberOperator (shift=1)
    test_mul(n_op1, "NumberOperator(shift=1)")
    test_mul5(n_op1, "NumberOperator(shift=1)")

    # DestroyPowerOperator for k = 2, 3
    for k in [2, 3]
        if k < N
            ak = DestroyPowerOperator(N, k)
            adk = adjoint(ak)
            test_mul(ak, "DestroyPowerOperator(k=$k)")
            test_mul5(ak, "DestroyPowerOperator(k=$k)")
            test_mul(adk, "adjoint(DestroyPowerOperator(k=$k))")
            test_mul5(adk, "adjoint(DestroyPowerOperator(k=$k))")
        end
    end
end

# ── Algebraic simplification tests ──────────────────────────────────────────
println("\n--- Algebraic Simplifications ---")

a = DestroyOperator(10)
ad = adjoint(a)

# a' * a → NumberOperator
result = ad * a
@assert result isa NumberOperator "Expected NumberOperator, got $(typeof(result))"
@assert result.shift == 0
println("  ✓ a' * a → NumberOperator(shift=0)")

# a * a' → NumberOperator(shift=1)
result2 = a * ad
@assert result2 isa NumberOperator "Expected NumberOperator, got $(typeof(result2))"
@assert result2.shift == 1
println("  ✓ a * a' → NumberOperator(shift=1)")

# a^k → DestroyPowerOperator
for k in [2, 3, 5]
    pk = a^k
    @assert pk isa DestroyPowerOperator "Expected DestroyPowerOperator, got $(typeof(pk))"
    @assert pk.k == k
end
println("  ✓ a^k → DestroyPowerOperator")

# a^1 → DestroyOperator (identity case)
@assert (a^1) === a
println("  ✓ a^1 → DestroyOperator (identity)")

# (a')^k → AdjointOperator{DestroyPowerOperator}
for k in [2, 3]
    pk = ad^k
    @assert pk isa AdjointOperator "Expected AdjointOperator, got $(typeof(pk))"
    @assert pk.L isa DestroyPowerOperator
    @assert pk.L.k == k
end
println("  ✓ (a')^k → adjoint(DestroyPowerOperator)")

# a * a → DestroyPowerOperator(N, 2)
@assert (a * a) isa DestroyPowerOperator
@assert (a * a).k == 2
println("  ✓ a * a → DestroyPowerOperator(k=2)")

# a * a^2 → a^3,  a^2 * a → a^3,  a^2 * a^3 → a^5
@assert (a * DestroyPowerOperator(10, 2)).k == 3
@assert (DestroyPowerOperator(10, 2) * a).k == 3
@assert (DestroyPowerOperator(10, 2) * DestroyPowerOperator(10, 3)).k == 5
println("  ✓ destroy × destroy compositions add powers correctly")

# a' * a' → adjoint(DestroyPowerOperator(N, 2))
adag2 = ad * ad
@assert adag2 isa AdjointOperator
@assert adag2.L isa DestroyPowerOperator
@assert adag2.L.k == 2
println("  ✓ a' * a' → adjoint(DestroyPowerOperator(k=2))")

# Verify compositions give correct numerical results
v = randn(ComplexF64, 10)
w1 = similar(v)
w2 = similar(v)
w_tmp = similar(v)
mul!(w_tmp, a, v)       # a * v
mul!(w1, a, w_tmp)      # a * (a * v)
mul!(w2, a * a, v)      # (a^2) * v  via DestroyPowerOperator
@assert isapprox(w1, w2; atol = 1.0e-12)
println("  ✓ (a*a)*v == a*(a*v) numerically")

# adjoint(adjoint(a)) == a  (double adjoint unwrap)
@assert adjoint(ad) === a
println("  ✓ adjoint(adjoint(a)) === a")

# NumberOperator is Hermitian: adjoint returns itself
n_op = NumberOperator(10)
@assert adjoint(n_op) === n_op
println("  ✓ adjoint(NumberOperator) === NumberOperator")

# Type-preserving constructors and promotion
a32 = DestroyOperator{Float32}(10)
a64 = DestroyOperator{Float64}(10)
@assert DestroyOperator(10) isa DestroyOperator{Float64}
@assert NumberOperator(10) isa NumberOperator{Float64}
@assert DestroyPowerOperator(10, 2) isa DestroyPowerOperator{Float64}
@assert (adjoint(a32) * a64) isa NumberOperator{Float64}
@assert (a32 * adjoint(a64)) isa NumberOperator{Float64}
@assert (a32 * a64) isa DestroyPowerOperator{Float64}
@assert (a32^2) isa DestroyPowerOperator{Float32}
@assert (adjoint(a32) * adjoint(a64)).L isa DestroyPowerOperator{Float64}
println("  ✓ operator constructors default to Float64 and mixed products promote T")

# ── concretize tests ─────────────────────────────────────────────────────────
println("\n--- concretize tests ---")

N = 8
a = DestroyOperator(N)
ad = adjoint(a)

# DestroyOperator concretize matches QuantumToolbox's spdiagm pattern
A_sp = concretize(a)
@assert A_sp ≈ spdiagm(1 => ComplexF64[sqrt(i) for i in 1:(N - 1)])
println("  ✓ concretize(DestroyOperator) correct")

# Create operator
A_sp_dag = concretize(ad)
@assert A_sp_dag ≈ spdiagm(-1 => ComplexF64[sqrt(i) for i in 1:(N - 1)])
println("  ✓ concretize(adjoint(DestroyOperator)) correct")

# NumberOperator
@assert concretize(NumberOperator(N)) ≈ spdiagm(0 => ComplexF64[i - 1 for i in 1:N])
@assert concretize(NumberOperator(N; shift = 1)) ≈ spdiagm(0 => ComplexF64[i for i in 1:N])
println("  ✓ concretize(NumberOperator) correct for shift=0 and shift=1")

# Verify a' * a concretize == num concretize
@assert concretize(ad * a) ≈ concretize(NumberOperator(N))
println("  ✓ concretize(a' * a) == concretize(NumberOperator)")

a32_sp = concretize(DestroyOperator{Float32}(N))
@assert eltype(a32_sp) == Float32
@assert concretize(DestroyPowerOperator{ComplexF32}(N, 2)) |> eltype == ComplexF32
println("  ✓ concretize preserves operator precision")

# ── Type stability tests ────────────────────────────────────────────────────
println("\n--- Type Stability ---")

N = 10
a = DestroyOperator(N)
v = randn(ComplexF64, N)
w = similar(v)

@inferred mul!(w, a, v)
println("  ✓ mul!(w, DestroyOperator, v) is type-stable")

@inferred mul!(w, adjoint(a), v)
println("  ✓ mul!(w, adjoint(DestroyOperator), v) is type-stable")

@inferred mul!(w, NumberOperator(N), v)
println("  ✓ mul!(w, NumberOperator, v) is type-stable")

@inferred mul!(w, DestroyPowerOperator(N, 2), v)
println("  ✓ mul!(w, DestroyPowerOperator, v) is type-stable")

@inferred mul!(w, adjoint(DestroyPowerOperator(N, 2)), v)
println("  ✓ mul!(w, adjoint(DestroyPowerOperator), v) is type-stable")

# ── SciMLOperators integration: AddedOperator and ScaledOperator ─────────────
println("\n--- SciMLOperators Integration ---")

N = 10
a = DestroyOperator(N)
ad = adjoint(a)

# a + a' → AddedOperator (no custom intercept, this is expected)
sum_op = a + ad
@assert sum_op isa SciMLOperators.AddedOperator
println("  ✓ a + a' → AddedOperator (expected)")

# verify a + a' mul! is correct
v = randn(ComplexF64, N)
w = similar(v)
sum_op_cached = cache_operator(sum_op, v)
mul!(w, sum_op_cached, v)
w_ref = concretize(a) * v + concretize(ad) * v
@assert isapprox(w, w_ref; atol = 1.0e-12)
println("  ✓ (a + a') * v gives correct result")

# scaled operator: 2.0 * a
scaled = 2.0 * a
v = randn(ComplexF64, N)
w = similar(v)
scaled_cached = cache_operator(scaled, v)
mul!(w, scaled_cached, v)
w_ref = 2.0 * concretize(a) * v
@assert isapprox(w, w_ref; atol = 1.0e-12)
println("  ✓ (2.0 * a) * v gives correct result")

# ── Float32 type preservation ────────────────────────────────────────────────
println("\n--- Float32 Type Preservation ---")

N = 10
a = DestroyOperator(N)
v32 = randn(ComplexF32, N)
w32 = similar(v32)

mul!(w32, a, v32)
@assert eltype(w32) == ComplexF32
println("  ✓ mul!(w::ComplexF32, DestroyOperator, v::ComplexF32) preserves ComplexF32")

mul!(w32, adjoint(a), v32)
@assert eltype(w32) == ComplexF32
println("  ✓ mul!(w::ComplexF32, adjoint(DestroyOperator), v::ComplexF32) preserves ComplexF32")

mul!(w32, NumberOperator(N), v32)
@assert eltype(w32) == ComplexF32
println("  ✓ mul!(w::ComplexF32, NumberOperator, v::ComplexF32) preserves ComplexF32")

mul!(w32, DestroyPowerOperator(N, 2), v32)
@assert eltype(w32) == ComplexF32
println("  ✓ mul!(w::ComplexF32, DestroyPowerOperator, v::ComplexF32) preserves ComplexF32")

# Verify Float32 results are numerically close to Float64 results
v64 = ComplexF64.(v32)
w64 = similar(v64)
mul!(w32, a, v32)
mul!(w64, a, v64)
@assert isapprox(ComplexF64.(w32), w64; atol = 1.0e-6)
println("  ✓ Float32 and Float64 results agree (within Float32 precision)")

# ── NormalOrderedOperator tests ──────────────────────────────────────────────
println("\n--- NormalOrderedOperator: mul! correctness ---")

for N in [5, 10, 20, 50]
    for (k, n) in [(2, 1), (1, 2), (2, 3), (3, 2), (2, 2), (3, 3), (1, 3), (3, 1)]
        if max(k, n) < N
            L = NormalOrderedOperator(N, k, n)
            test_mul(L, "NormalOrderedOperator(N=$N, k=$k, n=$n)")
            test_mul5(L, "NormalOrderedOperator(N=$N, k=$k, n=$n)")
        end
    end
end

# NormalOrderedOperator adjoint: adjoint((â†)^k â^n) = (â†)^n â^k
println("\n--- NormalOrderedOperator: adjoint ---")

N = 10
for (k, n) in [(2, 3), (3, 2), (1, 3), (2, 2)]
    L = NormalOrderedOperator(N, k, n)
    Ladj = adjoint(L)
    @assert Ladj isa NormalOrderedOperator
    @assert Ladj.k == n && Ladj.n == k "adjoint should swap k and n"
    # Verify adjoint concretize matches transpose of original
    @assert concretize(Ladj) ≈ transpose(concretize(L))
end
println("  ✓ adjoint(NormalOrderedOperator) swaps k,n and concretize matches transpose")

# Self-adjoint when k == n
L_sym = NormalOrderedOperator(10, 2, 2)
L_sym_adj = adjoint(L_sym)
@assert concretize(L_sym) ≈ concretize(L_sym_adj)
println("  ✓ NormalOrderedOperator(k=n) is self-adjoint")

# Adjoint mul! correctness
for (k, n) in [(2, 3), (3, 1), (1, 2)]
    L = NormalOrderedOperator(N, k, n)
    Ladj = adjoint(L)
    test_mul(Ladj, "adjoint(NormalOrderedOperator(k=$k, n=$n))")
    test_mul5(Ladj, "adjoint(NormalOrderedOperator(k=$k, n=$n))")
end

# NormalOrderedOperator concretize
println("\n--- NormalOrderedOperator: concretize ---")

N = 10
# Compare against explicit matrix product: concretize((â†)^k) * concretize(â^n)
for (k, n) in [(2, 1), (1, 2), (2, 3), (3, 2), (2, 2)]
    L = NormalOrderedOperator(N, k, n)
    a = DestroyOperator(N)
    ad = adjoint(a)
    expected = concretize(ad^k) * concretize(a^n)
    @assert concretize(L) ≈ expected "concretize mismatch for k=$k, n=$n"
end
println("  ✓ concretize(NormalOrderedOperator) matches (â†)^k * â^n sparse product")

# NormalOrderedOperator algebraic simplification
println("\n--- NormalOrderedOperator: algebraic simplifications ---")

a = DestroyOperator(10)
ad = adjoint(a)

# (a')^k * a^n → NormalOrderedOperator
result = (ad^2) * (a^3)
@assert result isa NormalOrderedOperator "Expected NormalOrderedOperator, got $(typeof(result))"
@assert result.k == 2 && result.n == 3
println("  ✓ (a')^2 * a^3 → NormalOrderedOperator(k=2, n=3)")

result = (ad^3) * (a^2)
@assert result isa NormalOrderedOperator
@assert result.k == 3 && result.n == 2
println("  ✓ (a')^3 * a^2 → NormalOrderedOperator(k=3, n=2)")

# (a')^k * a → NormalOrderedOperator(k, 1)
result = (ad^2) * a
@assert result isa NormalOrderedOperator
@assert result.k == 2 && result.n == 1
println("  ✓ (a')^2 * a → NormalOrderedOperator(k=2, n=1)")

# a' * a^n → NormalOrderedOperator(1, n)
result = ad * (a^3)
@assert result isa NormalOrderedOperator
@assert result.k == 1 && result.n == 3
println("  ✓ a' * a^3 → NormalOrderedOperator(k=1, n=3)")

# Verify numerical correctness of algebraic simplification
v = randn(ComplexF64, 10)
w1 = similar(v)
w2 = similar(v)
w_tmp = similar(v)

# (a')^2 * a^3: compare direct NormalOrderedOperator vs sequential application
L = (ad^2) * (a^3)
mul!(w1, L, v)
# Sequential: first a^3 * v, then (a')^2 * result
a3 = a^3
ad2 = ad^2
mul!(w_tmp, a3, v)
mul!(w2, ad2, w_tmp)
@assert isapprox(w1, w2; atol = 1.0e-12)
println("  ✓ NormalOrderedOperator mul! matches sequential (â†)^2 â^3 application")

# a' * a still gives NumberOperator (preserved)
result = ad * a
@assert result isa NumberOperator
println("  ✓ a' * a still gives NumberOperator (not NormalOrderedOperator)")

# Type promotion
a32 = DestroyOperator{Float32}(10)
a64 = DestroyOperator{Float64}(10)
result = (adjoint(a32)^2) * (a64^3)
@assert result isa NormalOrderedOperator{Float64}
println("  ✓ NormalOrderedOperator type promotion works correctly")

# NormalOrderedOperator type stability
println("\n--- NormalOrderedOperator: type stability ---")

N = 10
L = NormalOrderedOperator(N, 2, 3)
v = randn(ComplexF64, N)
w = similar(v)

@inferred mul!(w, L, v)
println("  ✓ mul!(w, NormalOrderedOperator, v) is type-stable")

# NormalOrderedOperator Float32
println("\n--- NormalOrderedOperator: Float32 ---")

v32 = randn(ComplexF32, N)
w32 = similar(v32)
L32 = NormalOrderedOperator{Float32}(N, 2, 3)
mul!(w32, L32, v32)
@assert eltype(w32) == ComplexF32
println("  ✓ NormalOrderedOperator{Float32} preserves ComplexF32")

println("\n" * "="^60)
println("Matrix (AbstractVecOrMat) mul! Tests")
println("="^60)

# ── Test that mul!(W_mat, op, V_mat) matches column-wise mul! ───────────────
println("\n--- Matrix mul! for all BosonicOperators ---")

N = 10
ncols = 4
V = randn(ComplexF64, N, ncols)
W_mat = similar(V)
W_col = similar(V)

for (label, op) in [
        ("DestroyOperator", DestroyOperator(N)),
        ("adjoint(DestroyOperator)", adjoint(DestroyOperator(N))),
        ("NumberOperator(shift=0)", NumberOperator(N)),
        ("NumberOperator(shift=1)", NumberOperator(N; shift = 1)),
        ("DestroyPowerOperator(k=2)", DestroyPowerOperator(N, 2)),
        ("adjoint(DestroyPowerOperator(k=2))", adjoint(DestroyPowerOperator(N, 2))),
        ("NormalOrderedOperator(k=2,n=3)", NormalOrderedOperator(N, 2, 3)),
    ]
    # 3-arg
    mul!(W_mat, op, V)
    for c in 1:ncols
        mul!(view(W_col, :, c), op, view(V, :, c))
    end
    @assert isapprox(W_mat, W_col; atol = 1.0e-12) "$label: matrix 3-arg mul! mismatch"
    println("  ✓ $label — 3-arg matrix mul! matches column-wise")

    # 5-arg
    W_mat2 = randn(ComplexF64, N, ncols)
    W_col2 = copy(W_mat2)
    α, β = randn(ComplexF64), randn(ComplexF64)
    mul!(W_mat2, op, V, α, β)
    for c in 1:ncols
        mul!(view(W_col2, :, c), op, view(V, :, c), α, β)
    end
    @assert isapprox(W_mat2, W_col2; atol = 1.0e-12) "$label: matrix 5-arg mul! mismatch"
    println("  ✓ $label — 5-arg matrix mul! matches column-wise")
end

println("\n" * "="^60)
println("KroneckerOperator Tests")
println("="^60)

using SciMLOperators: IdentityOperator, MatrixOperator

# ── Helper: reference Kronecker product via concretize ───────────────────────
function kron_ref(L::KroneckerOperator)
    return concretize(L)
end

# ── 2-mode basic test ────────────────────────────────────────────────────────
println("\n--- KroneckerOperator: 2-mode basic ---")

NA, NB = 4, 5
a = DestroyOperator(NA)
b = DestroyOperator(NB)

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
K_ref = kron(concretize(a), concretize(b))

v = randn(ComplexF64, NA * NB)
w = similar(v)
mul!(w, K, v)
@assert isapprox(w, K_ref * v; atol = 1.0e-12)
println("  ✓ KroneckerOperator(a, b) matches kron(concretize(a), concretize(b))")

# With cache
K_cached = cache_operator(K, v)
w2 = similar(v)
mul!(w2, K_cached, v)
@assert isapprox(w2, K_ref * v; atol = 1.0e-12)
println("  ✓ cached KroneckerOperator gives same result")

# ── Single operator per mode ─────────────────────────────────────────────────
println("\n--- KroneckerOperator: single active operator ---")

# Only mode 1 active
K1 = KroneckerOperator((NA, NB), 1 => a)
K1_ref = kron(concretize(a), sparse(1.0I, NB, NB))
mul!(w, K1, v)
@assert isapprox(w, K1_ref * v; atol = 1.0e-12)
println("  ✓ A ⊗ I matches kron(concretize(A), I)")

# Only mode 2 active (innermost → no permutation fast path)
K2 = KroneckerOperator((NA, NB), 2 => b)
K2_ref = kron(sparse(1.0I, NA, NA), concretize(b))
mul!(w, K2, v)
@assert isapprox(w, K2_ref * v; atol = 1.0e-12)
println("  ✓ I ⊗ B matches kron(I, concretize(B)) — fast path (no permutation)")

# ── 4-mode with identities ──────────────────────────────────────────────────
println("\n--- KroneckerOperator: 4-mode with identities ---")

dims4 = (3, 4, 5, 6)
a4 = DestroyOperator(dims4[1])
b4 = DestroyOperator(dims4[4])
total4 = prod(dims4)

K4 = KroneckerOperator(dims4, 1 => a4, 4 => b4)
K4_ref = kron(
    concretize(a4), sparse(1.0I, dims4[2], dims4[2]),
    sparse(1.0I, dims4[3], dims4[3]), concretize(b4)
)

v4 = randn(ComplexF64, total4)
w4 = similar(v4)
mul!(w4, K4, v4)
@assert isapprox(w4, K4_ref * v4; atol = 1.0e-10)
println("  ✓ A ⊗ I ⊗ I ⊗ B (4-mode) matches materialized kron")

# Cached version
K4c = cache_operator(K4, v4)
w4c = similar(v4)
mul!(w4c, K4c, v4)
@assert isapprox(w4c, K4_ref * v4; atol = 1.0e-10)
println("  ✓ cached A ⊗ I ⊗ I ⊗ B gives same result")

# Different active positions: modes 1, 3
K4b = KroneckerOperator(dims4, 1 => a4, 3 => DestroyOperator(dims4[3]))
K4b_ref = kron(
    concretize(a4), sparse(1.0I, dims4[2], dims4[2]),
    concretize(DestroyOperator(dims4[3])), sparse(1.0I, dims4[4], dims4[4])
)
mul!(w4, K4b, v4)
@assert isapprox(w4, K4b_ref * v4; atol = 1.0e-10)
println("  ✓ A ⊗ I ⊗ B ⊗ I (4-mode) matches materialized kron")

# All 4 active
a1 = DestroyOperator(dims4[1])
a2 = NumberOperator(dims4[2])
a3 = DestroyOperator(dims4[3])
a4_all = DestroyOperator(dims4[4])
K4all = KroneckerOperator(dims4, 1 => a1, 2 => a2, 3 => a3, 4 => a4_all)
K4all_ref = kron(concretize(a1), concretize(a2), concretize(a3), concretize(a4_all))
mul!(w4, K4all, v4)
@assert isapprox(w4, K4all_ref * v4; atol = 1.0e-10)
println("  ✓ A ⊗ B ⊗ C ⊗ D (all active, 4-mode) matches materialized kron")

# ── 5-arg mul! ───────────────────────────────────────────────────────────────
println("\n--- KroneckerOperator: 5-arg mul! ---")

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
K_ref_mat = kron(concretize(a), concretize(b))
K_cached = cache_operator(K, v)
v = randn(ComplexF64, NA * NB)

# β = 0 case
w = randn(ComplexF64, NA * NB)
α = randn(ComplexF64)
mul!(w, K_cached, v, α, 0.0)
@assert isapprox(w, α * K_ref_mat * v; atol = 1.0e-12)
println("  ✓ mul!(w, K, v, α, 0) matches α * K_ref * v")

# β ≠ 0 case
w = randn(ComplexF64, NA * NB)
w_orig = copy(w)
α, β = randn(ComplexF64), randn(ComplexF64)
mul!(w, K_cached, v, α, β)
@assert isapprox(w, α * K_ref_mat * v + β * w_orig; atol = 1.0e-12)
println("  ✓ mul!(w, K, v, α, β) matches α * K_ref * v + β * w_orig")

# 5-arg with 4-mode
K4c = cache_operator(K4, v4)
w4 = randn(ComplexF64, total4)
w4_orig = copy(w4)
α, β = randn(ComplexF64), randn(ComplexF64)
mul!(w4, K4c, v4, α, β)
@assert isapprox(w4, α * K4_ref * v4 + β * w4_orig; atol = 1.0e-10)
println("  ✓ 5-arg mul! with 4-mode KroneckerOperator matches reference")

# ── All-identity (M = 0) ────────────────────────────────────────────────────
println("\n--- KroneckerOperator: all-identity (M=0) ---")

K_id = KroneckerOperator((IdentityOperator(3), IdentityOperator(4)))
v_id = randn(ComplexF64, 12)
w_id = similar(v_id)
mul!(w_id, K_id, v_id)
@assert isapprox(w_id, v_id; atol = 1.0e-15)
println("  ✓ all-identity KroneckerOperator copies input")

# 5-arg all-identity
w_id2 = randn(ComplexF64, 12)
w_id2_orig = copy(w_id2)
α, β = randn(ComplexF64), randn(ComplexF64)
mul!(w_id2, K_id, v_id, α, β)
@assert isapprox(w_id2, α * v_id + β * w_id2_orig; atol = 1.0e-12)
println("  ✓ 5-arg all-identity KroneckerOperator: α*v + β*w")

# ── Adjoint ──────────────────────────────────────────────────────────────────
println("\n--- KroneckerOperator: adjoint ---")

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
Kadj = adjoint(K)
K_ref_mat = kron(concretize(a), concretize(b))
Kadj_ref = K_ref_mat'

v = randn(ComplexF64, NA * NB)
w = similar(v)
Kadj_cached = cache_operator(Kadj, v)
mul!(w, Kadj_cached, v)
@assert isapprox(w, Kadj_ref * v; atol = 1.0e-12)
println("  ✓ adjoint(KroneckerOperator) matches adjoint of materialized kron")

# ── concretize ───────────────────────────────────────────────────────────────
println("\n--- KroneckerOperator: concretize ---")

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
K_mat = concretize(K)
K_ref_mat = kron(concretize(a), concretize(b))
@assert isapprox(K_mat, K_ref_mat; atol = 1.0e-12)
println("  ✓ concretize(KroneckerOperator(a,b)) matches kron")

K4_mat = concretize(K4)
@assert isapprox(K4_mat, K4_ref; atol = 1.0e-10)
println("  ✓ concretize(4-mode KroneckerOperator) matches reference")

# ── Full tuple constructor ───────────────────────────────────────────────────
println("\n--- KroneckerOperator: full tuple constructor ---")

K_full = KroneckerOperator((a, IdentityOperator(NB)))
K_sparse = KroneckerOperator((NA, NB), 1 => a)
v = randn(ComplexF64, NA * NB)
w_full = similar(v)
w_sparse = similar(v)
mul!(w_full, K_full, v)
mul!(w_sparse, K_sparse, v)
@assert isapprox(w_full, w_sparse; atol = 1.0e-15)
println("  ✓ full tuple constructor matches sparse constructor")

K_full3 = KroneckerOperator((a, IdentityOperator(NB), b))
K_sparse3 = KroneckerOperator((NA, NB, NB), 1 => a, 3 => b)
v3 = randn(ComplexF64, NA * NB * NB)
w3a = similar(v3)
w3b = similar(v3)
mul!(w3a, K_full3, v3)
mul!(w3b, K_sparse3, v3)
@assert isapprox(w3a, w3b; atol = 1.0e-12)
println("  ✓ full tuple constructor with 3 modes matches sparse constructor")

# ── Base.kron overload ───────────────────────────────────────────────────────
println("\n--- KroneckerOperator: kron overload ---")

# kron with KroneckerOperator flattens
K_base = KroneckerOperator((NA, NB), 1 => a, 2 => b)
K_flat = kron(K_base, DestroyOperator(3))
@assert K_flat isa KroneckerOperator
@assert length(K_flat.dims) == 3
println("  ✓ kron(KroneckerOperator, op) flattens correctly")

K_flat2 = kron(DestroyOperator(3), K_base)
@assert K_flat2 isa KroneckerOperator
@assert length(K_flat2.dims) == 3
println("  ✓ kron(op, KroneckerOperator) flattens correctly")

# kron(K, K) flattening
K_flat3 = kron(K_base, K_base)
@assert K_flat3 isa KroneckerOperator
@assert length(K_flat3.dims) == 4
println("  ✓ kron(KroneckerOperator, KroneckerOperator) flattens correctly")

# Numerical verification of flattened kron
v_flat = randn(ComplexF64, NA * NB * 3)
w_flat = similar(v_flat)
K_flat_cached = cache_operator(K_flat, v_flat)
mul!(w_flat, K_flat_cached, v_flat)
K_flat_ref = kron(concretize(a), concretize(b), concretize(DestroyOperator(3)))
@assert isapprox(w_flat, K_flat_ref * v_flat; atol = 1.0e-12)
println("  ✓ flattened kron(K, op) gives correct numerical result")

# ── Mixed operator types ────────────────────────────────────────────────────
println("\n--- KroneckerOperator: mixed operator types ---")

M_mat = randn(ComplexF64, NB, NB)
M_op = MatrixOperator(M_mat)
K_mixed = KroneckerOperator((NA, NB), 1 => a, 2 => M_op)
K_mixed_cached = cache_operator(K_mixed, v)
v = randn(ComplexF64, NA * NB)
w = similar(v)
mul!(w, K_mixed_cached, v)
K_mixed_ref = kron(concretize(a), M_mat)
@assert isapprox(w, K_mixed_ref * v; atol = 1.0e-12)
println("  ✓ BosonicOperator ⊗ MatrixOperator gives correct result")

# ── size and properties ──────────────────────────────────────────────────────
println("\n--- KroneckerOperator: size and properties ---")

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
@assert size(K) == (NA * NB, NA * NB)
@assert size(K, 1) == NA * NB
@assert islinear(K) == true
@assert has_adjoint(K) == true
println("  ✓ size, islinear, has_adjoint correct")

# ── Type stability ───────────────────────────────────────────────────────────
println("\n--- KroneckerOperator: type stability ---")

K = KroneckerOperator((NA, NB), 1 => a, 2 => b)
K_cached = cache_operator(K, randn(ComplexF64, NA * NB))
v = randn(ComplexF64, NA * NB)
w = similar(v)

@inferred mul!(w, K_cached, v)
println("  ✓ mul!(w, KroneckerOperator, v) is type-stable (cached)")

@inferred mul!(w, K_cached, v, 1.0, 0.0)
println("  ✓ mul!(w, KroneckerOperator, v, α, β) is type-stable (cached)")

println("\n" * "="^60)
println("All tests passed!")
println("="^60)
