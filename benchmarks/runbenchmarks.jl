using BenchmarkTools
using QuantumToolbox

BLAS.set_num_threads(1)

const SUITE = BenchmarkGroup()

include("correlations_and_spectrum.jl")
include("dynamical_fock_dimension.jl")
include("eigenvalues.jl")
include("steadystate.jl")
include("timeevolution.jl")

benchmark_correlations_and_spectrum!(SUITE)
benchmark_dfd!(SUITE)
benchmark_eigenvalues!(SUITE)
benchmark_steadystate!(SUITE)
benchmark_timeevolution!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE, verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
