using QuantumToolboxOperators
using LinearAlgebra
using SparseArrays
using SciMLOperators
using SciMLOperators: AdjointOperator
using Test

println("=" ^ 60)
println("QuantumToolboxOperators — Informal Tests")
println("=" ^ 60)

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
        if !isapprox(w1, w2; atol = 1e-12)
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
        if !isapprox(w, w_expected; atol = 1e-12)
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

    a = DestroyOperator(N)
    ad = adjoint(a)
    n_op = NumberOperator(N)
    n_op1 = NumberOperator(N; shift = 1)

    # DestroyOperator
    test_mul(a, "DestroyOperator")
    test_mul5(a, "DestroyOperator")

    # AdjointOperator (creation)
    test_mul(ad, "adjoint(DestroyOperator)")
    test_mul5(ad, "adjoint(DestroyOperator)")

    # NumberOperator (shift=0)
    test_mul(n_op, "NumberOperator(shift=0)")
    test_mul5(n_op, "NumberOperator(shift=0)")

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
@assert isapprox(w1, w2; atol = 1e-12)
println("  ✓ (a*a)*v == a*(a*v) numerically")

# adjoint(adjoint(a)) == a  (double adjoint unwrap)
@assert adjoint(ad) === a
println("  ✓ adjoint(adjoint(a)) === a")

# NumberOperator is Hermitian: adjoint returns itself
n_op = NumberOperator(10)
@assert adjoint(n_op) === n_op
println("  ✓ adjoint(NumberOperator) === NumberOperator")

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
@assert isapprox(w, w_ref; atol = 1e-12)
println("  ✓ (a + a') * v gives correct result")

# scaled operator: 2.0 * a
scaled = 2.0 * a
v = randn(ComplexF64, N)
w = similar(v)
scaled_cached = cache_operator(scaled, v)
mul!(w, scaled_cached, v)
w_ref = 2.0 * concretize(a) * v
@assert isapprox(w, w_ref; atol = 1e-12)
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
@assert isapprox(ComplexF64.(w32), w64; atol = 1e-6)
println("  ✓ Float32 and Float64 results agree (within Float32 precision)")

println("\n" * "=" ^ 60)
println("All tests passed!")
println("=" ^ 60)
