using BenchmarkTools
using QuantumToolbox

suite = BenchmarkGroup()

suite["steadystate"] = BenchmarkGroup(["steadystate"])


## steadystate ##

N = 50
Δ = 0.1
F = 2
γ = 1
a = destroy(N)
H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ) * a]

suite["steadystate"]["driven-dissipative harmonic oscillator"] = @benchmarkable steadystate($H, $c_ops)

## end ##


BenchmarkTools.tune!(suite)
results = run(suite, verbose = true)
display(median(results))

BenchmarkTools.save(joinpath(@__DIR__, "benchmarks_output.json"), median(results))