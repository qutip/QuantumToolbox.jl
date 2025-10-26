using Revise
using LinearAlgebra
using SparseArrays
using QuantumToolbox
# using GenericLinearAlgebra
using GenericSchur
using BenchmarkTools

# include("arnoldi_new.jl")

# %%

# T = ComplexF64
T = Complex{BigFloat}

# A = rand(T, 100, 100)
A = sprand(T, 100, 100, 0.01) + Diagonal(rand(T, 100))
x = rand(T, 100)

nev = 5
krylov_dim = 15
tol = 1e-8

# values0 = eigvals(Array(A); sortby=abs2)[end-nev+1:end] |> reverse
res0 = eigen(Array{ComplexF64}(A); sortby=abs2)
values0 = res0.values[end-nev+1:end]
vectors0 = res0.vectors[:, end-nev+1:end]

res = eigsolve(A; v0=x, eigvals=nev, krylovdim=krylov_dim, eigstol=tol);
res.iter

values0
res.values

idxs1 = sortperm(res.values, by=abs2)
idxs2 = sortperm(values0, by=abs2)

res.values[idxs1] - values0[idxs2]

# res.vectors[:, idxs1] - vectors0[:, idxs2]

values_2 = [v' * A * v for v in eachcol(res.vectors[:, idxs1])]

values_2 - values0[idxs2]

map(1:nev) do i
    v0 = vectors0[:, idxs2[i]] |> copy
    v = res.vectors[:, idxs1[i]] |> copy

    abs(dot(v0, v))
end

# %%

@benchmark eigsolve($A; v0=$x, eigvals=$nev, krylovdim=$krylov_dim, eigstol=$tol)
