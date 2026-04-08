using LinearAlgebra
using SparseArrays
using SciMLOperators
using QuantumToolboxOperators
# using CUDA
using Metal
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

# a_gpu = adapt(CuArray, a)
a_gpu = adapt(MtlArray, a)

ψ_gpu = adapt(MtlArray, ψ)
dψ_gpu = adapt(MtlArray, dψ)

mul!(dψ_gpu, a_gpu, ψ_gpu)

@be mul!(dψ_gpu, a_gpu, ψ_gpu)

@benchmark mul!($dψ_gpu, $a_gpu, $ψ_gpu)

# %%

# a_reactant = Reactant.to_rarray(a)
a_reactant = adapt(Reactant.ConcreteRArray, a)

ψ_reactant = Reactant.to_rarray(ψ)
dψ_reactant = Reactant.to_rarray(dψ)

mul_compiled! = @compile mul!(dψ_reactant, a_reactant, ψ_reactant)

mul_compiled!(dψ_reactant, a_reactant, ψ_reactant)
