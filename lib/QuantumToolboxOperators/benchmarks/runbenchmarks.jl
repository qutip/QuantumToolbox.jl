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

# %%

N = 1000000
T = ComplexF32
a = DestroyOperator{T, true}(N)
# a_sparse = spdiagm(1 => sqrt.(T.(1:N-1)))
a_sparse = concretize(a)

ψ = randn(T, N) |> normalize
dψ = similar(ψ)

mul!(dψ, a, ψ)
mul!(dψ, a_sparse, ψ)

# %%

@be mul!(dψ, a, ψ)
@be mul!(dψ, a_sparse, ψ)

@benchmark mul!(dψ, a, ψ)
@benchmark mul!(dψ, a_sparse, ψ)

# %%

Adapt.@adapt_structure DestroyOperator
Adapt.@adapt_structure NumberOperator
Adapt.@adapt_structure DestroyPowerOperator
Adapt.@adapt_structure SciMLOperators.ScaledOperator
Adapt.@adapt_structure SciMLOperators.AdjointOperator
Adapt.@adapt_structure SciMLOperators.ComposedOperator
Adapt.@adapt_structure SciMLOperators.AddedOperator

# Adapt.adapt_structure(A::SciMLOperators.ScaledOperator) = SciMLOperators.ScaledOperator(A.λ, Adapt.adapt_structure(A.A))
# Adapt.adapt_structure(A::SciMLOperators.AdjointOperator) = SciMLOperators.AdjointOperator(Adapt.adapt_structure(A.L))
# Adapt.adapt_structure(A::SciMLOperators.ComposedOperator) = SciMLOperators.ComposedOperator(Adapt.adapt_structure.(A.ops), Adapt.adapt_structure.(A.cache))

a_gpu = adapt(CuArray, a)
# a_gpu = adapt(MtlArray, a)

a_sparse_gpu = CUSPARSE.CuSparseMatrixCSR(a_sparse)

ψ_gpu = adapt(CuArray, ψ)
dψ_gpu = adapt(CuArray, dψ)
# ψ_gpu = adapt(MtlArray, ψ)
# dψ_gpu = adapt(MtlArray, dψ)

mul!(dψ_gpu, a_gpu, ψ_gpu)

@be mul!(dψ_gpu, a_gpu, ψ_gpu)
@be mul!(dψ_gpu, a_sparse_gpu, ψ_gpu)

@benchmark mul!($dψ_gpu, $a_gpu, $ψ_gpu)
@benchmark mul!($dψ_gpu, $a_sparse_gpu, $ψ_gpu)

# %%

a_reactant = Reactant.to_rarray(a)
# a_reactant = adapt(Reactant.ConcreteRArray, a)

ψ_reactant = Reactant.to_rarray(ψ)
dψ_reactant = Reactant.to_rarray(dψ)

mul_compiled! = @compile mul!(dψ_reactant, a_reactant, ψ_reactant)

mul_compiled!(dψ_reactant, a_reactant, ψ_reactant)

@be mul_compiled!(dψ_reactant, a_reactant, ψ_reactant)

# %%

Δ = 0.1f0
U = 0.2f0
F = 0.3f0

H = Δ * (a' * a) + U * (a'^2 * a^2) + F * (a + a')
H = cache_operator(H, ψ)
H_sparse = Δ * (a_sparse' * a_sparse) + U * (a_sparse'^2 * a_sparse^2) + F * (a_sparse + a_sparse')

H_gpu = adapt(CuArray, H)
H_gpu = cache_operator(H_gpu, ψ_gpu)
H_sparse_gpu = CUSPARSE.CuSparseMatrixCSR(H_sparse)

H_reactant = Reactant.to_rarray(H)
H_reactant = cache_operator(H_reactant, ψ_reactant)

# %%

mul!(dψ, H, ψ)
mul!(dψ, H_sparse, ψ)

@be mul!(dψ, H, ψ)
@be mul!(dψ, H_sparse, ψ)

mul!(dψ_gpu, H_gpu, ψ_gpu)
mul!(dψ_gpu, H_sparse_gpu, ψ_gpu)

@be mul!(dψ_gpu, H_gpu, ψ_gpu)
@be mul!(dψ_gpu, H_sparse_gpu, ψ_gpu)

@benchmark mul!($dψ_gpu, $H_gpu, $ψ_gpu)
@benchmark mul!($dψ_gpu, $H_sparse_gpu, $ψ_gpu)

mul_compiled! = @compile mul!(dψ_reactant, H_reactant, ψ_reactant)
mul_compiled!(dψ_reactant, H_reactant, ψ_reactant)

@be mul_compiled!(dψ_reactant, H_reactant, ψ_reactant)

@benchmark mul_compiled!($dψ_reactant, $H_reactant, $ψ_reactant)

# %%

Base.summarysize(H)
Base.summarysize(H_sparse)
