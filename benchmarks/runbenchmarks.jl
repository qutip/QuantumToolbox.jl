using BenchmarkTools
using QuantumToolbox

fib(n) = n <= 1 ?  1 : fib(n - 2) + fib(n - 1)

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


tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("benchmarks_output.json", mean(results))