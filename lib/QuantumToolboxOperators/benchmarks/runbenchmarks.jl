using LinearAlgebra
using SparseArrays
using SciMLOperators
using QuantumToolboxOperators
using CUDA
# using Metal
using Reactant
using Adapt
using Chairmarks
using BenchmarkTools
using Markdown

# %%

N = 1000000
T = ComplexF32
a = DestroyOperator{T}(N)
# a_sparse = spdiagm(1 => sqrt.(T.(1:N-1)))
a_sparse = concretize(a)

Base.summarysize(a_sparse) / Base.summarysize(a)

# %%

ψ = randn(T, N) |> normalize
dψ = similar(ψ)

ψ_gpu = CuArray(ψ)
# ψ_gpu = MtlArray(ψ)
dψ_gpu = similar(ψ_gpu)

ψ_reactant = Reactant.to_rarray(ψ)
dψ_reactant = similar(ψ_reactant)

# %%

a_sparse_gpu = CUSPARSE.CuSparseMatrixCSR(a_sparse)

# %%

# ------- CPU -------

mul!(dψ, a, ψ)
mul!(dψ, a_sparse, ψ)

@be mul!(dψ, a, ψ)
@be mul!(dψ, a_sparse, ψ)

# ------- GPU -------

mul!(dψ_gpu, a, ψ_gpu)
mul!(dψ_gpu, a_sparse_gpu, ψ_gpu)

@be mul!($dψ_gpu, $a, $ψ_gpu)
@be mul!($dψ_gpu, $a_sparse_gpu, $ψ_gpu)

# ------- Reactant -------

mul_compiled! = @compile mul!(dψ_reactant, a, ψ_reactant)
mul_compiled!(dψ_reactant, a, ψ_reactant)

@be mul_compiled!($dψ_reactant, $a, $ψ_reactant)

# %% -------------- Real Hamiltonian ---------------

Δ = 0.1f0
U = 0.2f0
F = 0.3f0

H = Δ * a' * a + U * (a'^2 * a^2) + F * (a + a')
H = cache_operator(H, ψ)
H_sparse = Δ * (a_sparse' * a_sparse) + U * (a_sparse'^2 * a_sparse^2) + F * (a_sparse + a_sparse')

H_sparse_gpu = CUSPARSE.CuSparseMatrixCSR(H_sparse)

# %%

mul!(dψ, H, ψ)
mul!(dψ, H_sparse, ψ)

@be mul!($dψ, $H, $ψ)
@be mul!($dψ, $H_sparse, $ψ)


mul!(dψ_gpu, H, ψ_gpu)
mul!(dψ_gpu, H_sparse_gpu, ψ_gpu)

@be mul!($dψ_gpu, $H, $ψ_gpu)
@be mul!($dψ_gpu, $H_sparse_gpu, $ψ_gpu)


mul_H_compiled! = @compile mul!(dψ_reactant, H, ψ_reactant)
mul_H_compiled!(dψ_reactant, H, ψ_reactant)

@be mul_H_compiled!($dψ_reactant, $H, $ψ_reactant)

# %%

memory_ratio = Base.summarysize(H_sparse) / Base.summarysize(H)

# %% -------------- Print all the results on a table ---------------

bench_a_cpu = (tmp = @be mul!(dψ, a, ψ); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
bench_a_sparse_cpu = (tmp = @be mul!(dψ, a_sparse, ψ); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6

bench_a_gpu = (tmp = @be mul!(dψ_gpu, a, ψ_gpu); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
bench_a_sparse_gpu = (tmp = @be mul!(dψ_gpu, a_sparse_gpu, ψ_gpu); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
# bench_a_sparse_gpu = missing

bench_a_reactant = (tmp = @be mul_compiled!($dψ_reactant, $a, $ψ_reactant); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6

bench_H_cpu = (tmp = @be mul!(dψ, H, ψ); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
bench_H_sparse_cpu = (tmp = @be mul!(dψ, H_sparse, ψ); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6

bench_H_gpu = (tmp = @be mul!(dψ_gpu, H, ψ_gpu); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
bench_H_sparse_gpu = (tmp = @be mul!(dψ_gpu, H_sparse_gpu, ψ_gpu); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6
# bench_H_sparse_gpu = missing

bench_H_reactant = (tmp = @be mul_H_compiled!($dψ_reactant, $H, $ψ_reactant); sum(x -> x.time, tmp.samples) / length(tmp.samples)) * 1.0e6

md"""
Memory ratio (sparse/dense): $(round(memory_ratio * 1e-3, digits=2)) k

| Operator | CPU (Lazy) | CPU (Sparse) | GPU (Lazy) | GPU (Sparse) | Reactant (Lazy) |
|:--------:|:----------:|:------------:|:----------:|:------------:|:----------------:|
| a        | $(round(bench_a_cpu, digits=2)) μs | $(round(bench_a_sparse_cpu, digits=2)) μs | $(round(bench_a_gpu, digits=2)) μs | $(bench_a_sparse_gpu === missing ? "N/A" : string(round(bench_a_sparse_gpu, digits=2)) * " μs") | $(round(bench_a_reactant, digits=2)) μs |
| H        | $(round(bench_H_cpu, digits=2)) μs | $(round(bench_H_sparse_cpu, digits=2)) μs | $(round(bench_H_gpu, digits=2)) μs | $(bench_H_sparse_gpu === missing ? "N/A" : string(round(bench_H_sparse_gpu, digits=2)) * " μs") | $(round(bench_H_reactant, digits=2)) μs |
"""
