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
println("All tests passed!")
println("="^60)
