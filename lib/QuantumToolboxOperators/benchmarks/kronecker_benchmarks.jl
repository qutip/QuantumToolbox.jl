using Revise
using LinearAlgebra
using SciMLOperators
using QuantumToolboxOperators
using Chairmarks

# %%

# ─── Setup ───────────────────────────────────────────────────────────────────

T = ComplexF32
dims = (30, 30, 30, 30)       # 4 subsystems, each dim 30 → total 810,000
N_total = prod(dims)

println("Total Hilbert space dimension: $N_total")
println()

v = randn(T, N_total) |> normalize
w = similar(v)

a1 = DestroyOperator{T}(dims[1])
a2 = DestroyOperator{T}(dims[2])
a3 = DestroyOperator{T}(dims[3])
a4 = DestroyOperator{T}(dims[4])

a1_sp = concretize(a1)
a2_sp = concretize(a2)
a3_sp = concretize(a3)
a4_sp = concretize(a4)

I1 = IdentityOperator(dims[1])
I2 = IdentityOperator(dims[2])
I3 = IdentityOperator(dims[3])
I4 = IdentityOperator(dims[4])

# ─── Case 1: A ⊗ B ⊗ C ⊗ D (all active) ────────────────────────────────────

println("="^60)
println("Case 1: A ⊗ B ⊗ C ⊗ D (all active, 4 modes)")
println("="^60)

K1 = LocalTensorProductOperator(dims, 1 => a1, 2 => a2, 3 => a3, 4 => a4)
K1c = cache_operator(K1, v)

S1 = kron(a1, a2, a3, a4)
S1c = cache_operator(S1, v)

println("  LocalTensorProductOperator memory: $(Base.summarysize(K1c)) bytes")
println("  SciML TensorProduct memory: $(Base.summarysize(S1c)) bytes")
println("  Memory ratio (SciML / LocalTensor): $(round(Base.summarysize(S1c) / Base.summarysize(K1c); digits = 1))x")
println()

# Verify correctness
w_k = similar(v); w_s = similar(v)
mul!(w_k, K1c, v)
mul!(w_s, S1c, v)
println("  Correctness check: $(isapprox(w_k, w_s; atol = 1000 * eps(real(T))) ? "PASS" : "FAIL")")

println("\n  LocalTensorProductOperator:")
display(@be mul!($w, $K1c, $v))
println("\n  SciML TensorProductOperator:")
display(@be mul!($w, $S1c, $v))
println()

# ─── Case 2: A ⊗ I ⊗ I ⊗ B ─────────────────────────────────────────────────

println("="^60)
println("Case 2: A ⊗ I ⊗ I ⊗ B (2 active, 2 identity)")
println("="^60)

K2 = LocalTensorProductOperator(dims, 1 => a1, 4 => a4)
K2c = cache_operator(K2, v)

S2 = reduce(kron, (a1, I2, I3, a4))
S2c = cache_operator(S2, v)

println("  LocalTensorProductOperator memory: $(Base.summarysize(K2c)) bytes")
println("  SciML TensorProduct memory: $(Base.summarysize(S2c)) bytes")
println("  Memory ratio (SciML / LocalTensor): $(round(Base.summarysize(S2c) / Base.summarysize(K2c); digits = 1))x")
println()

w_k = similar(v); w_s = similar(v)
mul!(w_k, K2c, v)
mul!(w_s, S2c, v)
println("  Correctness check: $(isapprox(w_k, w_s; atol = 1000 * eps(real(T))) ? "PASS" : "FAIL")")

println("\n  LocalTensorProductOperator:")
display(@be mul!($w, $K2c, $v))
println("\n  SciML TensorProductOperator:")
display(@be mul!($w, $S2c, $v))
println()

# ─── Case 3: A ⊗ I ⊗ B ⊗ I ─────────────────────────────────────────────────

println("="^60)
println("Case 3: A ⊗ I ⊗ B ⊗ I (alternating)")
println("="^60)

K3 = LocalTensorProductOperator(dims, 1 => a1, 3 => a3)
K3c = cache_operator(K3, v)

S3 = reduce(kron, (a1, I2, a3, I4))
S3c = cache_operator(S3, v)

println("  LocalTensorProductOperator memory: $(Base.summarysize(K3c)) bytes")
println("  SciML TensorProduct memory: $(Base.summarysize(S3c)) bytes")
println("  Memory ratio (SciML / LocalTensor): $(round(Base.summarysize(S3c) / Base.summarysize(K3c); digits = 1))x")
println()

w_k = similar(v); w_s = similar(v)
mul!(w_k, K3c, v)
mul!(w_s, S3c, v)
println("  Correctness check: $(isapprox(w_k, w_s; atol = 1000 * eps(real(T))) ? "PASS" : "FAIL")")

println("\n  LocalTensorProductOperator:")
display(@be mul!($w, $K3c, $v))
println("\n  SciML TensorProductOperator:")
display(@be mul!($w, $S3c, $v))
println()

# ─── Case 4: Single operator A ⊗ I ⊗ I ⊗ I ─────────────────────────────────

println("="^60)
println("Case 4: A ⊗ I ⊗ I ⊗ I (single active, fast path)")
println("="^60)

K4 = LocalTensorProductOperator(dims, 4 => a4)  # physics mode 4 → julia dim 1 (fast path)
K4c = cache_operator(K4, v)

S4 = reduce(kron, (I1, I2, I3, a4)) # kron(I1, I2, I3, a4)
S4c = cache_operator(S4, v)

println("  LocalTensorProductOperator memory: $(Base.summarysize(K4c)) bytes")
println("  SciML TensorProduct memory: $(Base.summarysize(S4c)) bytes")
println("  Memory ratio (SciML / LocalTensor): $(round(Base.summarysize(S4c) / Base.summarysize(K4c); digits = 1))x")
println()

w_k = similar(v); w_s = similar(v)
mul!(w_k, K4c, v)
mul!(w_s, S4c, v)
println("  Correctness check: $(isapprox(w_k, w_s; atol = 1000 * eps(real(T))) ? "PASS" : "FAIL")")

println("\n  LocalTensorProductOperator:")
display(@be mul!($w, $K4c, $v))
println("\n  SciML TensorProductOperator:")
display(@be mul!($w, $S4c, $v))
println()
