using LinearAlgebra
using SparseArrays
using QuantumToolboxOperators
using Chairmarks

# %%

N = 1000000
T = ComplexF32
a = DestroyOperator{T}(N)
a_sparse = spdiagm(1 => sqrt.(T.(1:N-1)))

ψ = randn(T, N) |> normalize
dψ = similar(ψ)

mul!(dψ, a, ψ)
mul!(dψ, a_sparse, ψ)

# %%

@be mul!(dψ, a, ψ)
@be mul!(dψ, a_sparse, ψ)
