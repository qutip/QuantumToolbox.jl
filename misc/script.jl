using LinearAlgebra
using QuantumToolbox

LinearAlgebra.ishermitian(L::QuantumToolbox.AddedOperator) = all(ishermitian, L.ops)

N = 300
Δ = 0.5
F = 1.0
γ = 0.1

a = destroy(N)

H = Δ * a' * a + F * (a + a')

c_ops =  (sqrt(γ) * a, )

L = liouvillian(H, c_ops; matrix_form = Val(true))

ρ = rand_dm(N)

L_cached = cache_operator(L, ρ)

res1 = L_cached.data * ρ.data

L_vec = liouvillian(H, c_ops; matrix_form = Val(false))

res2 = vec2mat(L_vec.data * vec(ρ.data))

res1 - res2

norm(res1 - res1')
norm(res2 - res2')

Base.summarysize(L_cached)
Base.summarysize(L_vec)

Base.summarysize(L_vec) / Base.summarysize(L_cached)

# %%

using Chairmarks

dρ = randn(ComplexF64, N, N)
mul!(dρ, L_cached.data, ρ.data)
mul!(vec(dρ), L_vec.data, vec(ρ.data))

@be mul!($dρ, $L_cached.data, $ρ.data)
@be mul!($(vec(dρ)), $L_vec.data, $(vec(ρ.data)))
