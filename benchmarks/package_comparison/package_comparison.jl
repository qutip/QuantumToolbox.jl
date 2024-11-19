#=
# Performance Comparison of Quantum Simulation Packages: Julia vs. Python

Here we compare the performance of [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) with other quantum simulation packages:
- [`QuTiP`](https://github.com/qutip/qutip) (Python)
- [`dynamiqs`](https://github.com/dynamiqs/dynamiqs) (Python - JAX)
- [`QuantumOptics.jl`](https://github.com/qojulia/QuantumOptics.jl) (Julia)

To allow reproducibility, this page is generated with [`Literate.jl`](https://github.com/fredrikekre/Literate.jl) based on [this envitonment](https://github.com/qutip/QuantumToolbox.jl/blob/main/benchmarks/package_comparison). Moreover, to keep the code clean, we use the [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl) package to call Python code from Julia. We tested that the overhead of calling Python code from Julia is negligible for the purpose of this benchmark.

## Importing the Required Packages
=#

import QuantumToolbox
import QuantumOptics
using CairoMakie
using PythonCall
using BenchmarkTools

np = pyimport("numpy")
qutip = pyimport("qutip")
jax = pyimport("jax")
jnp = jax.numpy
dynamiqs = pyimport("dynamiqs")

dynamiqs.set_device("cpu")

#=
## Master Equation simulation

Parameters:
=#

N = 50
Δ = 0.1
F = 2
γ = 1
nth = 0.8

#=
### QuantumToolbox.jl
=#

a = QuantumToolbox.destroy(N)
H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumToolbox.fock(N, 0)

QuantumToolbox.mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false)).states[2] # Warm-up

mesolve_quantumtoolbox = @benchmark QuantumToolbox.mesolve($H, $ψ0, $tlist, $c_ops, progress_bar = Val(false)).states[2]

#=
### QuTiP
=#

a = qutip.destroy(N)
H = Δ * a.dag() * a + F * (a + a.dag())
c_ops = pylist([np.sqrt(γ * (1 + nth)) * a, np.sqrt(γ * nth) * a.dag()])

tlist = np.linspace(0, 10, 100)
ψ0 = qutip.fock(N, 0)

qutip.mesolve(H, ψ0, tlist, c_ops).states[1] # Warm-up

mesolve_qutip = @benchmark qutip.mesolve($H, $ψ0, $tlist, $c_ops).states[1]

#=
### dynamiqs
=#

a = dynamiqs.destroy(N)
H = Δ * jnp.matmul(dynamiqs.dag(a), a) + F * (a + dynamiqs.dag(a))
c_ops = [jnp.sqrt(γ * (1 + nth)) * a, jnp.sqrt(γ * nth) * dynamiqs.dag(a)]

tlist = jnp.linspace(0, 10, 100)
ψ0 = dynamiqs.fock(N, 0)

dynamiqs.mesolve(H, c_ops, ψ0, tlist, options = dynamiqs.Options(progress_meter = nothing)).states # Warm-up

mesolve_dynamiqs =
    @benchmark dynamiqs.mesolve($H, $c_ops, $ψ0, $tlist, options = dynamiqs.Options(progress_meter = nothing)).states

#=
### QuantumOptics.jl
=#

bas = QuantumOptics.FockBasis(N)
a = QuantumOptics.destroy(bas)

H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumOptics.fockstate(bas, 0)

QuantumOptics.timeevolution.master(tlist, ψ0, H, c_ops)[2][2]

mesolve_quantumoptics = @benchmark QuantumOptics.timeevolution.master($tlist, $ψ0, $H, $c_ops)

#=
## Monte Carlo quantum trajectories simulation

Parameters:
=#

N = 50
Δ = 0.1
F = 2
γ = 1
nth = 0.8
ntraj = 100

#=
### QuantumToolbox.jl
=#

a = QuantumToolbox.destroy(N)
H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumToolbox.fock(N, 0)

QuantumToolbox.mcsolve(H, ψ0, tlist, c_ops, progress_bar = Val(false), ntraj = ntraj).states[2] # Warm-up

mcsolve_quantumtoolbox =
    @benchmark QuantumToolbox.mcsolve($H, $ψ0, $tlist, $c_ops, progress_bar = Val(false), ntraj = ntraj).states[2]

#=
### QuTiP
=#

a = qutip.destroy(N)
H = Δ * a.dag() * a + F * (a + a.dag())
c_ops = pylist([np.sqrt(γ * (1 + nth)) * a, np.sqrt(γ * nth) * a.dag()])

tlist = np.linspace(0, 10, 100)
ψ0 = qutip.fock(N, 0)

qutip.mcsolve(
    H,
    ψ0,
    tlist,
    c_ops,
    ntraj = ntraj,
    options = pydict(Dict("progress_bar" => false, "map" => "parallel", "num_cpus" => Threads.nthreads())),
).states[1] # Warm-up

mcsolve_qutip = @benchmark qutip.mcsolve(
    $H,
    $ψ0,
    $tlist,
    $c_ops,
    ntraj = ntraj,
    options = pydict(Dict("progress_bar" => false, "map" => "parallel", "num_cpus" => Threads.nthreads())),
).states[1]

#=
### dynamiqs (not yet implemented)

=#

#=
### QuantumOptics.jl
=#

bas = QuantumOptics.FockBasis(N)
a = QuantumOptics.destroy(bas)

H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumOptics.fockstate(bas, 0)

function quantumoptics_mcwf(tlist, ψ0, H, c_ops, ntraj)
    Threads.@threads for i in 1:ntraj
        QuantumOptics.timeevolution.mcwf(tlist, ψ0, H, c_ops, display_beforeevent = true, display_afterevent = true)[2][2]
    end
end

quantumoptics_mcwf(tlist, ψ0, H, c_ops, ntraj) # Warm-up

mesolve_quantumoptics = @benchmark quantumoptics_mcwf($tlist, $ψ0, $H, $c_ops, ntraj)

#=
## Plotting the Results
=#

mesolve_times = [
    1e-6 * sum(m.times) / length(m.times) for
    m in [mesolve_quantumtoolbox, mesolve_qutip, mesolve_dynamiqs, mesolve_quantumoptics]
]
mcsolve_times =
    [1e-6 * sum(m.times) / length(m.times) for m in [mcsolve_quantumtoolbox, mcsolve_qutip, mesolve_quantumoptics]]

##

fig = Figure(size = (500, 300))
ax = Axis(
    fig[1, 1],
    xticks = (1:2, ["mesolve", "mcsolve"]),
    ylabel = "Time (ms)",
    title = "Performance Comparison with Other Packages",
)

colors = Makie.wong_colors()

barplot!(
    ax,
    ones(length(mesolve_times)),
    mesolve_times,
    dodge = 1:length(mesolve_times),
    color = colors[1:length(mesolve_times)],
)

barplot!(ax, 2 * ones(length(mcsolve_times)), mcsolve_times, dodge = 1:length(mcsolve_times), color = colors[[1, 2, 4]])

ylims!(ax, 0, nothing)

# Legend
labels = ["QuantumToolbox.jl", "QuTiP", "dynamiqs", "QuantumOptics.jl"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

axislegend(ax, elements, labels, position = :lt)

save("package_comparison.png", fig);

## fig

#=
![](package_comparison.png)

---

## System Information
=#

using InteractiveUtils

versioninfo()

#=
---
=#

QuantumToolbox.about()

#=
---
=#

qutip.about()
