using BenchmarkTools
using LinearAlgebra
using SparseArrays
using QuantumToolbox
using SciMLBase: EnsembleSerial, EnsembleThreads
using ForwardDiff
using Zygote
using Mooncake
using Enzyme: Enzyme, Const, Active, Duplicated
using SciMLSensitivity: BacksolveAdjoint, EnzymeVJP, MooncakeVJP

BLAS.set_num_threads(1)

const SUITE = BenchmarkGroup()

include("correlations_and_spectrum.jl")
include("dynamical_fock_dimension.jl")
include("dynamical_shifted_fock.jl")
include("eigenvalues.jl")
include("steadystate.jl")
include("timeevolution.jl")
include("autodiff.jl")

benchmark_correlations_and_spectrum!(SUITE)
benchmark_dfd!(SUITE)
benchmark_dsf!(SUITE)
benchmark_eigenvalues!(SUITE)
benchmark_steadystate!(SUITE)
benchmark_timeevolution!(SUITE)
benchmark_autodiff!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE, verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
