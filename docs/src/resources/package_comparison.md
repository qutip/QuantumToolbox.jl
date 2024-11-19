```@meta
EditURL = "../../../benchmarks/package_comparison/package_comparison.jl"
```

# Performance Comparison of Quantum Simulation Packages: Julia vs. Python

Here we compare the performance of [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) with other quantum simulation packages:
- [`QuTiP`](https://github.com/qutip/qutip) (Python)
- [`dynamiqs`](https://github.com/dynamiqs/dynamiqs) (Python - JAX)
- [`QuantumOptics.jl`](https://github.com/qojulia/QuantumOptics.jl) (Julia)

To allow reproducibility, this page is generated with [`Literate.jl`](https://github.com/fredrikekre/Literate.jl) based on this source file: [`package_comparison.jl`](<unknown>/benchmarks/package_comparison/package_comparison.jl). Moreover, to keep the code clean, we use the [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl) package to call Python code from Julia. We tested that the overhead of calling Python code from Julia is negligible for the purpose of this benchmark.

## Importing the Required Packages

````julia
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
````

````
Python: None
````

## Master Equation simulation

Parameters:

````julia
N = 50
Δ = 0.1
F = 2
γ = 1
nth = 0.8
````

````
0.8
````

### QuantumToolbox.jl

````julia
a = QuantumToolbox.destroy(N)
H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumToolbox.fock(N, 0)

QuantumToolbox.mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false)).states[2] # Warm-up

mesolve_quantumtoolbox = @benchmark QuantumToolbox.mesolve($H, $ψ0, $tlist, $c_ops, progress_bar = Val(false)).states[2]
````

````
BenchmarkTools.Trial: 100 samples with 1 evaluation.
 Range (min … max):  49.161 ms … 56.433 ms  ┊ GC (min … max): 0.00% … 11.11%
 Time  (median):     49.788 ms              ┊ GC (median):    0.72%
 Time  (mean ± σ):   50.093 ms ±  1.048 ms  ┊ GC (mean ± σ):  0.69% ±  1.59%

     █▇▄▂                                                      
  ▆▆▆████▄▇▇▃▄▄▄▄▆▃▄▁▁▁▁▁▃▃▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃ ▃
  49.2 ms         Histogram: frequency by time        55.8 ms <

 Memory estimate: 7.42 MiB, allocs estimate: 802.
````

### QuTiP

````julia
a = qutip.destroy(N)
H = Δ * a.dag() * a + F * (a + a.dag())
c_ops = pylist([np.sqrt(γ * (1 + nth)) * a, np.sqrt(γ * nth) * a.dag()])

tlist = np.linspace(0, 10, 100)
ψ0 = qutip.fock(N, 0)

qutip.mesolve(H, ψ0, tlist, c_ops).states[1] # Warm-up

mesolve_qutip = @benchmark qutip.mesolve($H, $ψ0, $tlist, $c_ops).states[1]
````

````
BenchmarkTools.Trial: 22 samples with 1 evaluation.
 Range (min … max):  229.655 ms … 235.674 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     230.936 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   231.229 ms ±   1.442 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

      ▃ █      ▃█                                                
  ▇▁▁▇█▁█▁▇▇▇▁▇██▁▁▇▁▁▁▁▇▁▁▇▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  230 ms           Histogram: frequency by time          236 ms <

 Memory estimate: 464 bytes, allocs estimate: 20.
````

### dynamiqs

````julia
a = dynamiqs.destroy(N)
H = Δ * jnp.matmul(dynamiqs.dag(a), a) + F * (a + dynamiqs.dag(a))
c_ops = [jnp.sqrt(γ * (1 + nth)) * a, jnp.sqrt(γ * nth) * dynamiqs.dag(a)]

tlist = jnp.linspace(0, 10, 100)
ψ0 = dynamiqs.fock(N, 0)

dynamiqs.mesolve(H, c_ops, ψ0, tlist, options = dynamiqs.Options(progress_meter = nothing)).states # Warm-up

mesolve_dynamiqs =
    @benchmark dynamiqs.mesolve($H, $c_ops, $ψ0, $tlist, options = dynamiqs.Options(progress_meter = nothing)).states
````

````
BenchmarkTools.Trial: 14 samples with 1 evaluation.
 Range (min … max):  156.506 ms … 384.347 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     380.617 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   364.569 ms ±  59.936 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

                                                             ▄█  
  ▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄██ ▁
  157 ms           Histogram: frequency by time          384 ms <

 Memory estimate: 1.45 KiB, allocs estimate: 66.
````

### QuantumOptics.jl

````julia
bas = QuantumOptics.FockBasis(N)
a = QuantumOptics.destroy(bas)

H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumOptics.fockstate(bas, 0)

QuantumOptics.timeevolution.master(tlist, ψ0, H, c_ops)[2][2]

mesolve_quantumoptics = @benchmark QuantumOptics.timeevolution.master($tlist, $ψ0, $H, $c_ops)
````

````
BenchmarkTools.Trial: 56 samples with 1 evaluation.
 Range (min … max):  86.955 ms … 92.450 ms  ┊ GC (min … max): 0.00% … 0.99%
 Time  (median):     89.251 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   89.410 ms ±  1.343 ms  ┊ GC (mean ± σ):  0.14% ± 0.41%

    ▃      ▃   ▃ ▃▃█  ▃  █ ▃ █ ▃    ▃ ▃█▃ ▃                    
  ▇▁█▁▁▁▇▇▁█▇▇▁█▁███▇▇█▁▇█▇█▇█▇█▁▁▁▁█▁███▁█▇▁▁▇▇▁▇▇▁▁▁▁▁▇▁▇▁▇ ▁
  87 ms           Histogram: frequency by time        92.3 ms <

 Memory estimate: 4.78 MiB, allocs estimate: 643.
````

## Monte Carlo quantum trajectories simulation

Parameters:

````julia
N = 50
Δ = 0.1
F = 2
γ = 1
nth = 0.8
ntraj = 100
````

````
100
````

### QuantumToolbox.jl

````julia
a = QuantumToolbox.destroy(N)
H = Δ * a' * a + F * (a + a')
c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a']

tlist = range(0, 10, 100)
ψ0 = QuantumToolbox.fock(N, 0)

QuantumToolbox.mcsolve(H, ψ0, tlist, c_ops, progress_bar = Val(false), ntraj = ntraj).states[2] # Warm-up

mcsolve_quantumtoolbox =
    @benchmark QuantumToolbox.mcsolve($H, $ψ0, $tlist, $c_ops, progress_bar = Val(false), ntraj = ntraj).states[2]
````

````
BenchmarkTools.Trial: 99 samples with 1 evaluation.
 Range (min … max):  41.823 ms … 60.681 ms  ┊ GC (min … max): 0.00% … 7.26%
 Time  (median):     50.583 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   50.869 ms ±  4.035 ms  ┊ GC (mean ± σ):  2.02% ± 2.88%

                    ▆▃▃▆ ▁  █▁▆ ▁▃▃ ▁  ▁ ▆  ▁                  
  ▄▄▁▁▄▁▄▄▇▄▇▁▇▁▁▇▄▁████▄█▇▇███▇███▄█▄▁█▇█▇▄█▇▄▁▁▁▄▁▁▁▄▁▄▁▁▇▄ ▁
  41.8 ms         Histogram: frequency by time        60.6 ms <

 Memory estimate: 13.87 MiB, allocs estimate: 45362.
````

### QuTiP

````julia
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
````

````
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.125 s …   1.196 s  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.177 s              ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.170 s ± 28.192 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █                              █         █        █     █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁█▁▁▁▁▁█ ▁
  1.12 s         Histogram: frequency by time         1.2 s <

 Memory estimate: 1.95 KiB, allocs estimate: 51.
````

### dynamiqs (not yet implemented)

### QuantumOptics.jl

````julia
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
````

````
BenchmarkTools.Trial: 79 samples with 1 evaluation.
 Range (min … max):  54.977 ms … 73.196 ms  ┊ GC (min … max): 0.00% … 7.16%
 Time  (median):     63.828 ms              ┊ GC (median):    9.17%
 Time  (mean ± σ):   63.843 ms ±  3.861 ms  ┊ GC (mean ± σ):  9.87% ± 3.45%

        ▁  ▁             ▁▁▄▁▄ █▁▁▄▄   ▁█▁ ▄▄▁▄ ▁   ▁      ▁   
  ▆▁▁▆▁▆█▁▁█▆▁▆▆▁▆▆▆▁▆▆▁▆█████▆█████▁▆▆███▁████▆█▆▆▆█▆▁▆▁▁▆█▆ ▁
  55 ms           Histogram: frequency by time        70.8 ms <

 Memory estimate: 74.99 MiB, allocs estimate: 233869.
````

## Plotting the Results

````julia
mesolve_times = [
    1e-6 * sum(m.times) / length(m.times) for
    m in [mesolve_quantumtoolbox, mesolve_qutip, mesolve_dynamiqs, mesolve_quantumoptics]
]
mcsolve_times =
    [1e-6 * sum(m.times) / length(m.times) for m in [mcsolve_quantumtoolbox, mcsolve_qutip, mesolve_quantumoptics]]

#

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
````

Legend

````julia
labels = ["QuantumToolbox.jl", "QuTiP", "dynamiqs", "QuantumOptics.jl"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

axislegend(ax, elements, labels, position = :lt)
````

````
Makie.Legend()
````

save("package_comparison.png", fig)

````julia
fig
````

```@raw html
<img width=500 height=300 style='object-fit: contain; height: auto;' src="data:image/png;base64, iVBORw0KGgoAAAANSUhEUgAAA+gAAAJYCAIAAAB+fFtyAAAABmJLR0QA/wD/AP+gvaeTAAAgAElEQVR4nOzdeUBN6f8H8Oe2uS1SWlUqRZLIFllCMiZZIlszMrINZgZf+24Y+9jGOtYvhmSLIYOsWSqpUDIkEmnf9/3e3x/nO8/vmXNvt6TF5f366/Sc557z3HPuOX3Oc55FIBaLCQAAAAAAfNoUGroAAAAAAABQNQTuAAAAAAByAIE7AAAAAIAcQOAOAAAAACAHELgDAAAAAMgBBO4AAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgCBOwAAAACAHEDgDgAAAAAgBxC4AwAAAADIAQTuAAAAAAByAIE7AAAAAIAcQOAOAAAAACAHELgDAAAAAMgBBO4AAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgCBOwAAAACAHEDgDgAAAAAgBxC4AwAAAADIAQTuAAAAAAByAIE7AAAAAIAcQOAOAAAAACAHELgDAAAAAMgBBO4AAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgCBOwAAAACAHEDgDgAAAAAgB5QaugAA9ef27dvR0dFSVwmFQgsLCysrK0NDw/opjFgsDgsLi4mJKSws1NDQGDZsmFAorJ9dg2yRkZFXrlx5+fJlamqqtra2lZVVq1at+vfvr6Oj09BF++RkZWXZ29tzy5aWlv7+/g1bngZRs4Pg7e2dl5fHLU+ZMkVRUbG2yvPu3bsLFy78/fffSUlJSkpKzZo169ix4+DBg/X19Wuwtbor50eKjo6+ffu21FWKioqmpqZWVlampqb1XGAfH5+cnBxueeLEiSoqKvW5d/giiAG+GOPHj6/yiujQocO9e/fquiSpqandu3dn9xsfH1/XO4UqBQUF8c4Lpamp+csvv+Tl5TV0GT8taWlp9BC1bt26oYvTMGp2EMzMzOinioqKaqUkL168GDRokNQfsIKCwvjx4xMSEj50m3VRzlpx+PDhKu/nBgYGhw8fFolE9VYqS0tLuvesrKx62y98OdBUBuBfnjx54ujouGvXrjrdy9KlS4ODg+t0F/Ch9u/f7+joWNl5yc3NXbFiRbt27dgoDUC26dOnt/xHeHh4ne7r5MmT7du3/+uvv6SuFYlER48etbGxCQgIaNhy1qeUlJQJEyaMGjWqoQsCUGsQuANIsWjRori4uLrb/s2bN+mys7PzTz/9pK6uXne7gyqdPXt26tSpFRUVbKKWlpaenh6bEhcXN2rUqPLy8votHcirpKSk1/8oKiqqux35+Ph8++23paWlNEVBQcHY2NjQ0FAgENDEnJwcFxeXu3fvNlQ5G4Svr++ZM2cauhQAtQOBO3yhpkyZkvSP+Pj4oKCg6dOn07UFBQVbt26tu72/f/+eW1BRUfHz89u5c6e2tnbd7Q5kS0tL+/777+mfAoFg0aJFL168yMzMTE1NTU5OnjFjBl17586dNWvWNEQxP0U6OjpZ/wgNDW3o4jSMBj8Ir169+v7778ViMfenoqLi6tWr09LS3r9/n5SUlJiYOGfOHJq5pKRk7NixmZmZ9V/OutO2bVt6P09MTHz06NGvv/7auHFjmmH16tUNWDyAWoTOqfCFUlNTY/uhmpiYdO/ePT8//9ixY1zKs2fPJD+VmZn55MmTjIyMJk2atGnTpnnz5pJ50tLSSkpKuGVDQ0MlJSWRSBQeHp6Wlubq6sql05pdVVVVVVVVqSWsqKiIiIhISEgoKSnR1dW1s7OTGtzn5ubm5uZyy9ra2lzN/cuXL6OiogYOHKiqqpqdnZ2fn89l0NXVFQqFxcXFwcHBGRkZzZo169q1q7KyMt1adHR0dHR0o0aN2rRpY2pqWtnRKygoePv27du3b0tLS42NjU1NTaX2e6PPJ4qKis2aNSOE5OXlhYaGZmdnW1hYWFtby+6P++LFi9jY2JKSkpYtW1pZWTVq1EhG5jdv3kRHR+fl5enp6XXs2LFJkyYyMvPs3r07KyuLWxYIBKdPnx45ciRda2BgsGPHDkLIzp07uZT9+/cvX75cstNbQUHB48ePU1JSFBUVDQ0NO3bsKLXMKSkpZWVl3LKxsbFAIEhPTw8NDS0uLra0tGzfvj3NKRKJwsLC4uPj9fX127Rpo6ury9tUcnIyrf43MTEhhGRlZYWFheXl5ZmZmXXs2FFBodLameqcxJKSEto0SE1NrWnTplz5w8LCbG1tzczMRCIR/XU1atSIDZXoV0hISIiNjc3MzDQ1NTU3N6+yj29KSkpUVFRGRoa6urqZmVnbtm3ZOmNK8kKrqKgIDQ1NSEgwNDS0sbGpzsNwVlZWQUEBt6ylpaWhoUFXscfWyMiIHsnCwkIa9XLHpMqDUE1FRUUhISGZmZnVuTpYS5cupQXQ1NT09fXt378/XWtoaLhlyxZ7e/vvvvuO++G9f/9+y5Yta9euretylpSUhIeHJyUlNWrUqEWLFm3btpWaR/ZvrDpFUlJSYu/nXH9cHR2dSZMmcSnR0dHl5eVKSv+KeTIyMt6+ffvu3btGjRqZmJiYm5vLPnf5+fmPHj1KS0tr2rSptbU1d0+rPq6DAbdM79UcsViclJQUFxeXmJjYtGlTrjAyOraKxeInT568efNGS0ure/fuqqqqIpEoMTGRW6uuri754xeLxZGRke/evSsrKzM2Nu7UqRN755f6ZWNjY2NjYzU0NFq0aGFqaio7P9SfBm1hD1Cv2M6ps2bNksxw+vRpmqFZs2bsqrS0tJEjR/Jita5du966dYu3kT59+tAMT58+jYyMbNeuHSHE0tIyJSWlc+fOnTt3phkUFRW5lLCwMLqFsrKydevW8UJPBQWFgQMHxsbG8nb3888/0zw7d+7Myclxc3Pj/uQ6vP744480w/nz548ePaqlpUVT9PT0vL29xWLx33//7eDgwO6xb9++z58/5+0uLi5u6tSpkv9RnJ2db9y4wctM1+ro6JSXly9fvlxNTY0mqqurb9u2raKiQvJE+Pn5mZubs9tXVFScPHlySkqKZOYnT5506tSJzaysrDxixIhXr15JZpbKyMiIfnb48OFS82RkZLDf+tmzZ+za9PT07777jndYVFVVZ8yYIdmflQ1fUlNT58yZw8YTtra24eHhYrH44sWLXCxOfwDff/89r7ubsbExzVBYWLhixQo2ijI2Nvb395f8LtU/ieyoHSNGjBCJREuXLuU+eOzYMbFY/ObNG5qhW7du7GfLysoOHTokGXh179798uXLUg/y06dPHR0deWG6np7enj17JPsX8i60M2fOsKGbQCCYMmVKdna21B1RmzZtoh+ZMWMGTS8pKWEfuh4/fkxXLV++nKYvX768soOwevXqzp07s9da69atO3fufOrUKS4Dr9Pntm3b2EteV1f3yJEjsgvPiYuLY+9L69evryzn1KlTaTYtLa2CgoK6K2dZWdmyZct4LQBNTU337dvHu+Sr/I1Vhu2camdnJ5khNTWV3Xt0dDRdFRQU9PXXX5N/U1FR8fLyYrNRWVlZXl5evMjVwsLi3LlzvJyVdU5lb9TW1tbp6elcekVFxYkTJ9q0acMrjI6OzvLly6V2b3348GHLli1pzqZNm3p7e7N9b8aOHcv7yKlTp3iPGU2aNJk3b15OTo7k9l++fDl27FjeM7+GhsbChQtTU1OlnwyoRwjc4QtSZeB++fJlmsHMzIymBwQEVFa5IhAIdu7cyW6EjSeuXr1KazEtLS3j4+OlboQQcvv2be7jL1++7NChQ2XZ1NTUfv/9d3Z37P+D7du39+rVi/4pGbi7uLhI3ezmzZul1oM2b948NzeX7uv9+/e8Nt8sRUXFa9eusWWjq3R0dKZMmSL1UytXrmQ/Ul5ePmHChMp2oaWlxXuW2LFjR2X1Uk2aNAkMDKzyV/Hy5Uv2UyEhIZXl3L9//9J/vHjxgqZfunRJxkB7pqamvHGK2MBdMnoghGhra2/ZskVqZfmwYcPYTbGB+9y5cyXzKyoq/vbbb+xHPugk8oKqJUuW0D+rDNw9PT0r2wsh5JdffuEd3rVr18qo0nN0dOSNiMJeaIsWLZL6qb59+8o++5GRkTRzp06daPrDhw/Z7ezbt4+uGjBgAE3nfmBSDwKt6+Whtws2IF64cKHUzCdOnJBdfrFYzPakNzAwyM/Pryznu3fv2IuFe6iri3K+e/eua9euUnMSQjw8PMrLy2nmKn9jlakycKfvUjhv3rzh0q9evSrjTZShoSHvl3bv3j32QuOZOnUqm1lq4O7j40MTjY2N3759S/MvXry4si0TQpycnMrKytjt+/j4SH3FsX37drrMBu4lJSUyhlOztbVNTk5mtx8ZGclWr/AYGxsnJSXJOClQDxC4wxekysCdbbs8YMAALrGgoICtkbWysho5cqSdnR1NUVRUZKNJNp7o0qULXa5OjbtIJOINR6ipqcmre1ZQUGDr/9jAnd0dkRa4U1LbHkhNX7NmDd2Xk5MTTVdSUmrZsqWNjQ37X8TIyIg9ntXZhVAoZP9N7tu3j12rqqpqbW3N1n3a2dkVFxdzmcPDw9mtOTg4uLu7s+2XzM3NS0pKZP8qvL29aX51dfXq/ZT+X3JyMu/diJGRES8ybt68OVe7yZFsMFDZ6ZC66v79+3RTbDxBc0o2vmLf53zQSWSDKjs7Ozawlh24Hz16lC1As2bNevfuzZ4agUDA1u5fuXKF960tLS159bXu7u7skWcvNBmHy9fXV/YZpMdQSUmJniY2DCKETJw4kUsXiUT0dGtra3MB6MfXuFdWfgMDg9LSUtnlZ5t1jRs3Tnbmnj170sxLliypo3IOHTqUrtLX13dzc+vTpw/7WmD//v00c5W/scpUGbjfv3+fZlBVVeVq+jMzM9lKCnV1dRsbm5YtW7Lf6Ntvv6UbKS4ubtGiBftlzczMeNN9nDx5kuaXDNwfPHhAry8tLa2nT5/SzAEBAex+dXV127Vrx/674R2ruLg4Xj0F3TL71o4N3NnOWqqqql999ZWrqyt7y2K/bEVFBVuXr6ys3KVLly5durB3YEdHRxknBeoBAnf4gsgO3G/dusXezubMmcOlr1q1iibOnDmTvufdvXs3TR8xYgTdDi+eaNSo0dChQ2fNmrV06VKah/4Pa9KkCVuGgwcP0g+qq6ufPHmS293bt2/ZgL5nz570I2zgzundu/f06dOnT5+emZkp/nfgLhAINm3alJmZmZeXt3TpUvZTurq6Fy5cKCoqio+PZ6uB6VfLyMig/2PMzc3fvXtH09n3vOyA9Oz2dXR0jh49mpmZWVBQ4O3tzf77oa+bc3Nz2arrrVu3ctFAYmJix44daTrXtkf87xj0zJkzXGJxcTFb18t7HyKJ/cdmY2MjO7Mkdl+WlpaPHj3i0u/evcu+peGCJA4buOvo6Pz5559FRUUJCQkDBw5kj1inTp0ePXpUVlYWFhbGBk/sN+JVBI4dOzYhIaGiouLhw4dWVlY0ffDgwTU7iZIT3LRp08bLy2vmzJlBQUHiygN3Npq8ePEilygSibjeApxp06bRU8aGCy4uLlybqPLy8n379rE/FfZtAHuhKSgoLF++PD4+vqys7P79++yz7syZM2WfQbbKmb74+vbbb9lv3bZtWy7977//pomjR4/mEmW8dqDt1gghvBcv7DlVU1PbvXt3RkZGfn7+tm3b2F3zGmVJonM/EUJWrVolO/N3331HM7MRWy2W886dOzSxZ8+etLVSSEgIfRIzMTGhQ8JX+RurjOzAPSYmhr1p0Ncpvr6+NNHNzY2+oHj8+DENTy0tLel2NmzYQPPb2NjQavu9e/fS9FatWtH8vMD97du3BgYG3J9CofDu3btsIdle72vXrqX/XH7//XeaPmnSJJrfy8uLprdu3To8PFwkEsXExLBdGggTuGdlZXF9BgghOjo6tIIpKSmJa8BJCBEIBBEREVx6VFQU3YijoyM9OMnJyez9gVdJD/UMgTt8QdjA3cjIqOc/unXrxguAhEJhXFwc9ynay0dPT483/wgbOtB7IptoYGAQExMjWZLKAnc2fDlw4AC7Kicnh60Vo/9c2cBdIBBcuHCBty82cPfw8GBXsZ0d2d09fvyYptN/eLdv3zb7B6/1BRsNsI3+2UPKaww6efJkuoq2ymX/E48aNYrNzw5gxw2gERERQVOGDBnCO1a0k5mhoaHk8WexTZZdXFxkZ+Zh2z4pKCiw7WfE/45gVFVVad0/G7jv3buX5n/69Cl7xF6+fElXsUESfaQU/ztw79GjB7v3t2/fspWXXO+IDz2JvKBq/vz5vAbKlcWsrVq1oumRkZE0XSQSOTs7c9fdlClTuMTjx4/TzObm5rw6ZvbJ2dnZmaazFxovOmc3+PXXX4tlOnv2LM1M3y9ZWFhw55S7RhQUFLjWwP/9739p5v/+97+yD4K42gHx9u3b2VW9e/emq86ePSu7/FxROcePH5edmT2YX331VV2Uc/jw4TSR661BsTcr+suv8jdWGfZ2oa6uTu/nPXr0sLS05PVHoi8Q1q9fTy+BBw8esBuk3fEFAgFtzMN+fd5DVL9+/egq+tqQDdzj4+NpX3NFRcXz58/zvoKbmxtXEisrK/ZnHxsbSzfSu3dvLpHr9c4lKikp0f9QYrG4tLTU2tqafoQG7uyLoy1btrC7Zg87/b/Avn5ka6PEYvGxY8foEQ4ODq7OCYI6glFl4AuVmJhI++BLmj9/Pne/Tk1NpeONdOrUKTk5mc3WqVMnGpxFRESwt07OkiVL2FhctuTk5FevXnHLurq6vIaJmpqaEydOpNXD9+7dY5uzc4YNG8a+pJbUt29f9k8dHZ309HTJVWxATwfA6du3r+TY9omJiQEBAX5+fpL5WcrKyrwJHdm2RtnZ2dwC29x8zJgxbP5evXotWrSIGxOD+//64sULurZTp068stna2nJTKSUnJ6ekpNBKL0lshS4dpaSa7t27R5cHDx7cunVrdm3v3r25WnNCSFFRUXh4uOS0rJUddmNjYzb2lXpGeGbPns3+aWpqOnz4cNrf+vnz5y1atPiYk2hgYLB27VoZjYNZ7Lt1BwcHV1fXQYMGubi4GBoa3rhxg5eZfSqbOXMmr6X7jBkzfvnlF65IwcHBFRUVkuP5sMEiqeTXVZn+/fsrKSlxA8gEBgYSQtLT07nIqV27dhYWFufPnxeJRKGhoc7OziEhIfSDlfUYqQF3d3de+ekxofefyrCHOi8vT3bm4uJiuix5GKtUnXLSC7NJkyZNmzZlf282NjZ0+cmTJ5Lb/6DfGKugoIA7d1I5OjqOHj2aW160aBGvR0RFRcWrV698fX3fvXvHpXDhESGkpKSEPpzb2tqy5SeELF26lFbq5+bm8pq4EEI8PT1pJ4pFixYNGzaMl+HPP//kpRQXFz979mzz5s1s8biFZ8+e0eVhw4axTxTKysrz5s1ja0M47E3SxsaGPRcmJiYqKircwP/0XLC/JV9fX1tbWzc3t4EDB3bv3t3T01N2rxWoNwjcAf6lUaNG69atozFQTEwMXeXv789r7Mh6/fq1ZCJvtBPZ6L8NQoi1tbVkRz1bW1upmau/OxmDJLLV+TIUFxffvHnT398/LCwsKiqqykCBo6mpWVnTTBZ9biGEsPWIhBCBQLB+/Xo2hT01q1atYqsSeV6/fi0jcKevkkklR1UGtsadHcaRsrW15QJ3buOSgXtlZ6Sap4MlWYB27drRwJ39fdbsJNra2lZ/PLgRI0bQ1+6FhYVnz549e/asQCDo2LGjq6vr8OHD2d+q7MOora1tbGzMnZrCwsL09HTJs8kbKLP6AykSQpo0aeLg4MC1h+aqYB88eMCt6tGjh5mZ2fnz5wkhDx8+dHZ2pqvs7Ow+dDRAGXjll9E7UJKent7z58+5Zal3IRZ71cjoUV2ZKsspEolobXFOTs6H3jA/6DdWTV5eXmwDLU56evqVK1du3Ljx+PHj6OhoduIqXiFFIhG3zLsjEUL69evHVrpLYt+5Xbx4cdWqVVIflmJiYi5fvhwQEBAZGRkXF0f3yMPWa7CPphypvYHZ081ricd68+aNSCRSUFBwdnZu0qRJTk4Ol/7s2bNnz56tW7dOW1t7wIABrq6u7u7u7JCp0CAQuMMXqnnz5mwjYIFAYG5u3qFDhwEDBrA1nW/fvq3mBqX+H/qg8IuNn6QGmmygQO+tNd5dDQQEBEycOJFtGMBRVlamA5N/DDruO6nGd/mgU9OjR4/K1rKxxfv377l/YNXc8sefslokWQC2C11SUhK3UOOT+EG/rmXLlhUUFOzcuZN9iSEWix89evTo0aM1a9Y4OzufOHGCix2rcxjpM1VOTo6Mx7CaGThwIBe4Z2VlPX/+nEbn3bt3pyNyhoSEFBQU0KcRGWFQPWvbti2t9mZrWKViM0gOQfjxkpOTq/na6uNvmCwNDY1u3bqxKfr6+nZ2dj179pR8M7l79+5FixbRke8p+uKF+qA7kmxPnz7ds2cP26idEFJUVDRnzhxuzCI2XerFyN7xJB8a2aFjpX5EhpKSkoSEhObNm2tpad26dWvcuHFsXw5CSFZW1qlTp06dOjVr1qwtW7ZMnDixOpuFOoLAHb5Q7u7uv/32W5XZ2EopBwcHXnsPlmQ7mQ+lqalJl6U240lISKDLNZ7kpcYSEhLc3NzoZE9WVlbOzs4dO3a0t7c/ePAg21W3xtj/RmlpaTKq68i/T80333zDe4vNkn1qunfvLhAI6Jvx8+fPjxgxQmrOHj160OrzP/74Y/To0Z/UKUtMTOTV37NF4vpr1sNJ5CgpKW3atGnRokXHjx8/d+5cYGAgr/nNzZs3hw0bFhQURCR++eybJU5dH8avv/6a9tUOCgqi7WF69Oihr6+voKAgEolCQkLCwsLot/h0AncnJyfal/Hq1asxMTFs1QMrICCA7X3o7Oxc64XR1dVVVFTkjpKBgcFPP/1UWc7KJp6rGUtLS8kmWFL5+fmxperevbujo2OHDh169OgxePBg9vgQiTtSDQrm5uZ24cIFbnnFihUeHh7skFPLli2jnVyFQqGzs3O3bt06duzYunVrtl6Jw36QN0Q9YZ7MWfr6+rTSfcmSJTKOOW0k06lTp2fPnt25c+fEiRMXL17ktQ7Nzs6eNGmSiYkJOygq1DME7gCysHfPZs2aLVu2rO72xQ6F8fz586KiIt59Njw8nC5Xc0LBWnTgwAEa8M2bN+/XX3+l45PIGJ/+g7D9AaKjo3kvfx88eMDVh6mrq3fs2JE9NQ4ODjNnzqzZTrW0tHr27ElHjtu4caPUwD01NTUkJIS+xeaKyj5asGeHooE+qftTFh4ezqtDDQsLo8tc+/t6OIksHR2dWbNmzZo1KyMjw8/P78KFC1euXKE1ssHBwc+ePWvbti3vMPJiguTkZBq4q6mpyRiEvsY6deqkr6/PBUP379/nBnHX19fnOhra2NhERUUlJyfTdkeampoy3uHUM1dX16ZNm3KTuZaXl69YsYIdNZwS/3vIcG7Mx1ovjIqKiqmpKfc+p7S0tE5vmDWzZcsWuuzn5zd48GD6p2RLOW6YSO6pPjo6mrc2MTGRtgtq06aN5FQYs2bN+u2335ycnAICAggh2dnZS5YsOXDgALe2rKyMTsbctGnTR48e0VsEOzwAWxi6zHbN57AXO2VlZUWb/o8ZM0Zqcz6p+vTp06dPn71794aEhFy4cMHX15dtdXPgwAEE7g3og7uAAHxRmjdvTttxhoSE8N4Cc736OB8f9+jp6dHxRrKzs/fs2cPb15EjR+ifsttW1gX2/9bgwYNpwJeZmcl2LvwYbCzOG9A9NDSUqxtzdHTkeuiyPUElC/D333/TU1NlG272/XVoaOi0adN476krKiq8vLxo1G5iYsJ1SnN0dKTtaq5du8b73/nnn3/SN84aGhrssH11YfPmzeyL/sjIyEuXLnHLioqK3BuJejiJhJCrV686/IP7Gevo6Hh5eZ0/f/7t27fsoDpcfR7bQ3fHjh28eXM2bNhAGxI4OjryZq2vFQKBgA6B6uvryz3b0NCcnrg//viDW+D6s9Z6MWpGXV2drUI+efLk999/z7tTZWdnu7m50SZAhJBFixbVoA9oddALMysri1eBXVhYSK9K3hhK9YZeAoqKiq6urjT95s2b9JmWUlVVpQM3xcbG8ir1J0+e7PgPqX2gV65cSQhhH4//+9//0rtEbGwsvc+0b9+efbA/d+6c5NbY2+O5c+fY1m6lpaXssJWU7JtkUFAQdy5oJ/uRI0dyl62jo2NhYaFAIHBwcFi/fn10dDT7Oo5XDQ/1DIE7gCwCgeD777/nlhMTE+fPn09vtfHx8YMHD3ZycnJycnJ2dpZsMVkD7HBpixYt2rZtW3Z2dkVFRVBQUJ8+fYqKirhVAwYMcHBw+PjdfRC2I9rBgwe5B5WYmJgRI0ZUOXBHNY0cOZLWpwYGBk6dOpV7Pf3s2TN2sMIhQ4YQQuzt7WkHx3Pnzp04cYJmCAoK6t69O3dqvLy82KESKtsv2xB23759vXr12rhx461bt65fv75z5842bdqwMwStXLmS+09sYGAwbdo0mu7q6vrXX3+VlJQUFhZ6e3uPGzeOrlqyZEmt97rjiYiIcHV1DQ0NTU9PP3fuHDfYH7dqwoQJ3IGth5NICLGwsAj5xy+//MJWHxYVFdGfMfmnv/Xw4cNpXWBycnLfvn0jIiJEIlFaWtrSpUvpkHYCgWDFihW1VUgeOkQMvZBp4E4bT9NV1R9Phg2OJZs31JbFixfTYbkJIQcOHOjcufOiRYt8fX1PnDjxn//8p0OHDuyoQV9//TVv0KpaLCd7RUyfPp22MCkoKJg+fbrTP7hBn+ofvQQqKip27NiRk5NTUVEREBBQ2Qyy06dPp8uenp63bt2qqKgoLS1ds2YNvSe0adOGHQWSx97eno5pIxKJuKFLyb8vxrCwsGvXrpWWlhYVFR08eJAdVYZq3cI+WW0AACAASURBVLo1HQK1rKysX79+169fz8zMDAwM7Nu3L+2gzBo3bhwdO3/NmjW0nr6iouL333/v2bMndy5oUysVFRXusr1///7ixYvpGEQikYjtosP+2KABNMgglAANosqZU6XKyspix1LQ09Pr3r07OwckIcTLy4vmZ4eXZifJY1U2jrtYLP7qq6/YLQsEAt4QGUKhkB0bno31pU42xI7j7uPjw65i62PS0tJoOvv2gM5sItkAWurYF9w86hyaqKOjwysVfV9MCFm4cCFN51W0E0J4gxj06tWLDnjMjsZICGnRokX37t3puPucI0eOSD0FPG/evKnmICF0Ek1OVlYWb3QOFRUVXnVsmzZt2EGa2Spndv5wtpUqnfGHc+zYMbqK/enKmImdUlVVff/+fc1OIm86eqnHjWZghzDnTeLbvHlzW1tbOk42p0OHDjT//fv3eZNxSrbHnTx5MrtrGRca+1qfN7B6ZdLT03k10HSGWra9E4dOXCX7IIjFYrYu3MDAYPjw4XQ6KraGlTdBxMKFC+kq3nwOlYmOjpbaPVFSu3btUlNTeR+v3XKyrefV1dXt7e3btGnDDqjCDlte5W+sMlXOnCrVqFGj2KOhoKAgdQwiOuuC5Mypqqqq7AWuoKDw119/0e1LzpwqFotfv37NDqv1xx9/cOm8dl9CoVBy2Jnu3bvTjbNzwcrAzpy6evVqmq6oqGhra9u5c2f24hIKhXQSkosXL7KXobq6urW1tbW1Ne8mLDkgPdQn1LgDVEFLS8vHx4cO0JGWlhYcHMyGtp6envv376+t3R0+fJj9tycWi9mhlw0MDE6ePFn9seFrkZeXF68JdWFhISHEwsLCw8ODJtJxi2tm8uTJP/zwA/vPg32V0bJlyzNnztCq6169em3evJn+R3zz5k1wcDAdT1pBQWHnzp28msXKmJubBwYGsiG1JAUFBU9PT17sq6Wlde7cOfa/dWlpKdtkpVOnTmfPnq3r6nZCCDurIqWpqXn69Gka39fPSSSEXLx4kQ3T4+Pjo6Ki2DbEurq67MxHPXv23LdvH9vxlK2YJ4SMGTNGajVkbdHR0WHbfKuoqNA/27VrxwY6tra2vOd2Gdg+rCkpKVxjodooL5+VlVVwcHCVL+Lc3Nzu378v2U+gdst54MAB+uRWUFAQGhr6/Plz2q/Xzs4uICCgHq4IqVauXMmeTZFIxN1gXVxc2JlW6SXQqFGjEydOsLF7UVERvcAVFBTWr1/PNrmRysLCgn0RsWDBAq793saNG9lsxcXF3FFavHgxvQe+fPmS/gvo2bNnZePcc+8hJc2bN2/s2LHcckVFRVRUVHh4OL24tLW1b9y4QbvvDxkyZM2aNfSzBQUFL168ePHiBXsTnj17tuSA9FCfELgDVK1///5RUVGenp5s3YxAIHBycjp37twff/xRi/+EjI2Nb9y4cejQId67V01NzSlTpjx//pyd47A+qampXb58mZ2BRUND48cffwwLC2NjvpUrV0od36CaFBQUdu/effPmzW7durH/X5s1a7Zhw4aoqCh2iENCyNy5c0NDQ9m25oSQRo0aeXp6hoSEyBjUQlKLFi2ePHmye/fuDh068Fapq6u7ubk9efLk2LFjkvVzPXv2fPr06YIFC3ijXBsbG//6668PHz6UMeJNLdqxY8e8efPYg9ajR4+HDx+y3e/q5yQSQpo1axYZGbl161bJh0xdXd3FixdHRETwfuFTpkx59uzZyJEjeXXt7dq18/PzO3nypIxZCGoFG7x27NiRnmglJSV21PkPmnfJ1dV1586dtTjiuwwmJiZBQUEnT57s27cv752PUCgcNGjQjRs3/vzzT3YMnzoqZ4sWLYKDg1etWsXrr2ltbb1jx47AwMD6OSBS2djY+Pv7szG6kZHRtm3bLl26xL7DmTp1Kl12cHCIjIycOXMm+05DRUWlf//+YWFhCxYsqM5+ly9fTo98cnLyL7/8QgiZMGHC/v372XuanZ3d5cuX161bR0dqz8jIoEMeEUKWLFny119/sQ8SJiYm27dvP3TokNT9CoXC48ePnz17lncl6urqLly4MDIysmfPnmz6kiVLAgMDR48ezWthKBAI+vbte+HCBbZ3LzSI/3WXBoDqKCoqevnyZW5urr6+vqmpae2OaCYpISEhISGhpKREV1e3VatWn0h/uOzs7JiYGEVFxTZt2tTpERCJRG/evElLS7OysmKnSZIqMzMzNja2pKTExMTE2Nj4I49VUlLSq1evMjIydHV1LS0tqxlniMXiV69epaamKigoGBoayh7OslaYmJjQEVeKioqEQmFBQcHz589LSkpatGghOZUjVW8nkStYQkJCcnKypqamiYlJlaeypKQkJiYmMzNTTU3N1NS0BvMEfYLy8vJyc3M1NTU1NDR4jYLqQnZ29suXL1NSUriforW1NW3rXJ/lFIlEMTExqampTZo0MTU1reuJJj5IYmJiXFyclpaWtbV19fvp5ubmvnjxQktLy9LSsgZTz0pVUVERFxeXlJTUvHnz6o89xR3Yli1bctMaJCYm0hdrEydOlBrHv3v3LiEhQUlJibusZJ9fsViclpaWkJBQVFRkZGRkbGzcUC9JgAeBOwCAvJIM3Bu2PABQR4qKiuisVSoqKryB3q9cuUJb7Gzbtu0///lPfZcP6ssnUYEHAAAAAJUpKirq0KED7Spw6dIlOiFgTk7OokWLaM66HnkWGhYCdwAAAIBPWtOmTQcOHEgnZ/Dw8Pj666+7deuWlJR04sSJlJQULn3UqFG8ZuvwmUFTGQAAeYWmMgBfjvz8fBcXFzoZqiQnJyc/P79q9mcAOYVRZQAAAAA+dRoaGnfv3vXz8xs0aJCJiQnXO7ZRo0YtW7Z0c3O7ffv2rVu3ELV/9lDjDgAgr5KTk+mQ0tWcfwcAPg/l5eXZ2dk6Ojr1ME4RfDoQuAMAAAAAyAE0lQEAAAAAkAMI3AEAAAAA5AACdwAAAAAAOYDAHQAAAABADiBwBwAAAACQA5g59fOBAaEAAAAAGladDtiIGncAAAAAADmAGvfPDQbmh1qXlZWVl5enra3duHHjhi4LAMBHycvLy8rKaty4sba2dkOXBT439dD2ATXuAAAAAAByAIE7AAAAAIAcQOAOAAAAACAHELgDAAAAAMgBBO4AAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgAzp35Z6mFOLwCoMcx8DAAAMqDGHQAAAABADqDG/UuEWj2ATw3ehgEAQJVQ4w4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgCBOwAAAACAHEDgDgAAAAAgBxC4AwAAAADIAQTuAAAAAAByAOO4AwAAQO3LDZ1LRGUNXQq+0tJSxZKSchWV3EaNGrosUjTuvF6gpN7QpYBPlwBz8Xw2uAlcZJ/Q6uQBgPqHaxM+P8nHhOKKkoYuhZwx8EhVEOo1dCmghurhTo6mMgAAAAAAcgBNZYBv9fWXK65G1/9+l33VarWLdf3vFwAAAEAuoMYdAAAAAEAOIHAHAAAAAJADCNwBAAAAAOQAAncAAAAAADmAwB0AAAAAQA4gcAcAAAAAkAMI3AEAAAAA5AACdwAAAAAAOYAJmEC+hYeH37p1KyEhISsrS09Pz8zMbODAgS1btmzocsmTxMTE3Nzc6uc3NjZu3Lhx3ZWHEFJeXv7q1StCiLm5uVAoJIQUFha+e/eOENKqVStFRcU63XttiY2NLS0t1dfXb9q0aUOXBQAAPgcI3EFe+fj4LF269M2bN5Kr7O3tN23a1KdPn/ovlTyaOXOmr69v9fP7+Ph4eHjUXXkIIcnJyW3atCGEBAcHOzg4EEIePnzo5ORECElLS9PV1a3TvddAVFTU8ePHCSE9evQYOnQolzhgwIDXr1+vX79+0aJFDVo6AAD4TCBwB/lTVFTk4eFx8eJF7k99ff127dppa2unpKQ8fvw4Pz8/NDS0b9++c+fO3bRpk0AgaNjSynD16tWoqChLS8vhw4c3dFngo7x48WLjxo2EkB9//JEG7gAAALULgTvImbKyMhcXl7t37xJC2rZtu3nz5q+++oq2nSguLj516tTixYuTkpK2bNmSm5u7f//+Bi2vLCdPnjx69OigQYMaNnDftWvXhg0beIl2dnaFhYWDBg367bffeKsMDQ3rq2gAAADw/xC4g5xZvnw5F7UPGzbs1KlTKioq7FqhUDh+/PhBgwYNGDDg8ePHBw4c6N27t6enZwMVVj5IDcS5NxUaGhroMFAd9vb2hw8fJoRwLXwAAADqAgJ3kCexsbGbN28mhJibm588eZIXtVO6uroXL15s1apVcXHx3LlzR44cyXVwBKgjZmZmXl5eDV0KAAD4zCFwB3mybdu2iooKQsiaNWsaNWokI6eJicmsWbM2btyYmpp69OjRqVOncukHDhwoKSnp16+fjY0N7yN3796NjIxs3ry5m5sbmy4Wiy9evOjn5xcfH5+RkaGvr29qajpq1ChnZ2feFgIDAx8/fty2bVsnJ6eMjIwDBw4EBARkZGSYmpp27979hx9+UFNT43Levn372bNnL168IITExcXt2rWLEDJhwgR1dfWIiIh79+5paGhIBoLFxcUHDx4khIwYMaJZs2Zc4rlz5xITEwcMGGBlZeXv73/06NHXr1+rqKh07tx51qxZLVq0IITExsbu2bMnNDS0uLjYwsLC3d191KhR1TvkVQsJCTl+/HhUVFRubq6urq6Dg8OECRPMzc2lZi4sLDx06NDNmzcTEhJUVFRatGgxZMiQ0aNHf1BXBG4jV69eTUlJadKkib29/ZQpUywtLWtcvOPHj2dnZxsYGEgelqdPn965c4cQMmbMGD09vcqKVFZWtm/fPkLI0KFDTU1Nq/9dAAAAqk8gFosbugxQO7jQR/YJrU6e1ddfrrgaXbtlq45lX7Va7WItO4+xsXFiYmL79u2fPHlSZaiXk5Njbm6enZ391VdfXbt2jUvU0tLKyck5cODA5MmTefl/+umn3bt3Ozs737hxgyamp6c7OTlFRUVJbt/JyenKlSvs88O8efO2bNkyadKkH374YejQoQkJCbzCBwcHN2/enBAyefLkQ4cO8TYYHx9vYmLy22+/zZ4929jY+P3797wM6enpXOx47969Xr16cYm9evUKDAz8448/goODf//9dza/lpbW7du3Y2Njx40bV1hYyK6aPn36nj17pB43joaGRkFBwZgxY06ePFlZnvLy8pUrV65fv14kErHpmpqav//++7fffsvLHxER4eHhwT2usAYOHHjkyBF9fX2a8v79e+5A0VFlAgICuFFlwsLCPD09eRtp3Ljx/v37eWPdVL94hw8fnjhxIiHk+vXr/fv3p+lFRUXt27d/9erV4MGD/fz8KjsOhJD8/HxuiMwrV664uLhwiS1btqz+qDLVuTYB5EvyMaG4oqShSyFnDDxSFYSV1hHAJ64e7uSocQe5ER0dnZiYSAhxd3evTgVtkyZNnJ2dfX19AwMDS0tLK2tXI9uYMWOioqIUFBTc3Nz69eunra2dlJR05syZhw8f3r59e9myZZs2beJ9JC0tbejQoQUFBb/88ku3bt3y8vKOHj3q5+eXkJAwZcqUq1evEkKGDx9ubm5+7ty5x48ft27dmmuF36RJkxqUkLN58+bIyMgxY8aMHj1aKBSeOHHC29s7Ozt76NChycnJTZs2Xbt2rZ2d3YsXL1atWpWSkrJ3796pU6fa2dnVeI+EkIkTJx47dowQYmtrO3ToUEtLy4iIiDNnziQlJY0dOzYnJ2f69Ok0899//92tW7eSkhJlZWUPDw8HB4eCgoKAgIDLly9fuXKla9euL168qE6LptGjR8fGxvbu3XvAgAGamppBQUG+vr55eXnffPNNRUXF2LFja1C8CRMmnD9/3s/Pb/LkyVFRURoaGlz60qVLX716paenx73oAAAAaFgI3EFuREf/7z2Avb19NT9ib2/v6+tbWFiYkJDANRr5IG/evLl16xYhZP369QsWLKDpc+fOHTJkyF9//XX9+nXJT128eNHIyCg0NJR26xwxYsTw4cP//PPPgICAsrIyZWXlQYMGDRo06NWrV48fP27ZsuWyZcs+tGw8kZGRq1atWrFiBfenq6trbm4u17zHyMgoPDyc64Hq5OTUsWPH7t27i8Xihw8ffkzg/ujRI27kcg8Pj8OHD9OYe9GiRYMHD3706NGqVau+++47dXV1Ln3x4sUlJSVaWlrnz5/v27cvlzh//vy9e/f+8MMPb9++3bVr17x586rcb2xs7OzZs7ds2cI9vM2YMeP27dsjRozIyspatmzZqFGjuCe0Dy3e/v37bW1t3759O3/+fO7FRVBQ0Pbt27lVBgYGNT5QAAAAtUWhoQsAUF3p6encAtumQjaaMy0trQZ7jIyM1NHR0dfXZ2uOCSECgYAbwPH169dSP/jzzz/zBmMZP348IaSkpIR7aVDrzM3NlyxZwqbQISbXrl3Ljhvj4OBgZGRECImPj/+YPf78889isbhp06YHDhxga8qbNWvGDcGZkpJCW+OEh4dz4+4vWLCARu2cadOmDRo0iBCyYcOG4uLiKvdraWnJG57fyclp9erVhJC4uDhvb+8aFI8QYmhoyP25b9++27dvFxcXT5w4USQSeXl5DRs27MMODQAAQN1A4A5yo6Tkf20lZXdLZdGcpaWlNdijm5tbenp6SkoK13yZlZWVRSppx6agoMCF6SwTExNuoY6avvXo0UNJ6V8v0GjvVckZZLlVvJbfH+rBgweEkMmTJ9OGJVTnzp25JvjBwcFsZmVl5R9++EFyU//5z38IIRkZGS9fvqxyvzNmzKDD9lOTJ0/mWv+Hh4fXoHic0aNHjx49WiwWT5o0ae7cudHR0ebm5jt27KiySAAAAPUDgTvIDR0dHW4hMzOzmh/JyMjgFnR1dT++ACUlJdHR0RcvXly4cOGqVasqy9a8efPqP1rUFhqmf9CqGsvNzeVegLRr105qBi49NjaW+5N7NWFqaiq1HT/dCM0vg9TmPY0aNWrfvj3d0YcWj9qzZ4+hoeGbN2/27NmjoKBw9OhRyWc2AACAhoLAHeQG18CDECI5LElluJwCgaDGwatYLD5z5szIkSMtLCzU1NSsra3d3Nx+/fXX/Pz8yj7yJUwsStsIVTbso5mZGSHkzZs3bP7KMuvr66uqqrL5ZeC2LMnCwoJu4UOLR+no6HDt2gkhEyZM6N27d5XlAQAAqDfonApyo0uXLkKhsLi4+ObNm3Rcdtlu3rxJCOnQoUN1BmyRbMSSn58/ePBgbhhvQoiRkVHbtm2trKy6dOmSl5c3c+bMD/wGH+vTGStQWVmZWygrK5OagWuHQwvM5a8ss1gs5nJW5wtWthGuJZWWllYNiseiHY4vX76clZWlra1dZZEAAADqB2rcQW4IhUKuubavr+/Tp0+l5omIiAgNDeWW/f39uYFoBg8eXJ3tS3bWnDdv3p07dwQCwcKFC+Pi4hISEq5du7Zr1y4vL68Giec+sjtpLeKqtwkhb9++lZohLi6OENKqVSvuT252pMoyJycnc91SaX4ZKquV52rZuT7BH1o8ys/P7+DBg4qKijo6OklJSfX/bAYAACADAneQJ7NnzyaEiESiOXPmSK5NT093c3Pr3bu3t7d3UVHR/PnzCSGqqqo//fQTL6fUqlY63CR14sQJQsi0adM2bNjAa6FRs2Fqqq+aJWwoampqXOsj+pjEw6XToXW4wP39+/dJSUmSmR8+fMgt8Ibiker+/fuSiZmZmU+ePKFb+NDicdLS0qZMmUIIWbBgATfyzPHjx7nBcAAAAD4FCNxBnnz99dfcYII3btzYsmULb62GhkavXr2Ki4s9PT27du3K1covXLiQHT6SG5CEq3BlXb9+nTekSWFhYV5eHiGE6/XII3sezY/BlTAlJUVybMTdu3fX0U5rwNnZmRBy5MgRyWcYf3//x48fE0K4uU4JIX379lVUVKyoqNi6davkpjZu3EgIMTExqU7gvnv3bskOBlu3bi0oKBAIBO7u7jUoHuf7779PSUmxsbH5+eef3d3dR4wYQQiZOnVq9TtDAwAA1CkE7iBnTpw4wc2GM2/ePE9Pz3fv3tFVQqHw2LFjo0aNIoRERUURQpycnHhzG1lZWRFCDh8+TAecIYS8efNmxowZvB2pqak1b96c/NNQniotLZ07d+7t27cJISUlJVxw/zF4YShXwoqKim3btrHp69evDwwM/Mh91aJVq1YpKysXFha6u7uzwXF4ePjkyZMJIa1atZo0aRKXaGVl5eXlRQjZsWPHkSNHaOby8vKZM2dywzKuWbOGtk2XISsra8iQITk5OTRl3759XOg/evRo+pT1QcUjhBw+fPjPP/9UVFQ8fPgwNyjQ7t27tbW1k5OTeb+NkJAQDw8PDw+PQ4cOVe9QAQAA1A50TgU506xZs3v37rm6ur569crb29vHx6dr1642NjZNmjRJTk4OCQlhx/h7//59TEyMtbU1TRkxYsSDBw+SkpJsbW2nTJliZGQUERFx6tSpwsLCadOm7d27l92Xh4fHpk2bzp496+LiMnDgQKFQGBMTc/r06fj4eEdHx3v37pWXl7u6uk6cOHHChAk1+C5NmzYlhNy7d8/FxUVXV3f79u06Ojrdu3dv1qxZUlLSkiVLQkNDnZyccnJybty4cefOnSFDhkRFRVVn6JV6YGFhMX/+/HXr1t2/f9/GxqZPnz4WFhaRkZF37twpLi5WVlbetm0bG4ivXLnSz88vNTV1woQJe/fu7datW35+/t27d1+9ekUI6dmz57hx46rcqbGxsY6OTkBAgJWVlZOTk46OTlBQENdIxtDQcM2aNTUrXlxc3KxZswghc+fO7dq1K5doYGCwdevWCRMmnDhxYvTo0W5ublx6fHz8qVOnCCFaWlps6A8AAFDXELiD/GnVqlVoaOi6det27txZXFz84MEDbrYdytLScuLEiYcOHYqJienWrZuPj4+rqyu3as6cOU+ePPH29k5OTuam2ySE6OjonD9/PiMjgxe4r1mz5sWLF35+fv7+/v7+/lyipqbmkSNHxo8fP2bMmNOnT9+/f19dXb1mgfv48eP3799fUFDAbXzDhg3c9k+dOjVq1KiUlJTz58+fP3+eyzxq1KjDhw9LHcW8oaxdu7Z169Y//vhjenq6r68vTW/ZsuWJEyfs7e3ZzCYmJpGRkd999921a9dCQkJCQkK4dIFAMHv27PXr1ysoVP0CUCgUXrp0acSIEaGhoVz0zHFwcPD19aUDhn5Q8UQi0fjx4/Py8tq0acMbnt/Ly8vHx+fatWvTpk1zdHTkHrQAAAAaiuDTGWAOPhI3CbzsE1qdPKuvv1xxtQE6QS77qtVqF+uq8zGys7OvXr169+7dxMTE/Px8IyMjc3NzBwcHFxcXBQWFjIyMkSNHBgQEDBo06NKlS+wHg4ODr1+//vr1a6FQaGdn5+Hh0bRp05SUlOfPn2tra/OC4zt37gQEBMTGxgqFwg4dOri7u3NtdcrLyy9cuJCcnNytW7cuXboQQl6/fh0fH6+pqdmpUydeUfPy8rhJPR0cHIRCIU1PSUnx9/fPy8szMjJydXWlMzfl5+d7e3v//fffmZmZRkZGQ4YM4Sb7DAkJKSoq6tixIx3g8vHjxzk5Oc2bN+c6gFKZmZmRkZGEkN69e/Ni4vDw8Ly8PDMzsxYtWlR2bO/du1dRUaGvr29jYyP7LCQkJJw/f/7p06d5eXm6uroODg5ubm7q6upSM4vF4uvXr9+8eTMxMVFZWdnCwmLo0KGSvQhKSkq49jOdOnXS1NQkhGRnZz958kRVVbVbt27l5eV//fXX7du309LSdHV1+/TpM2TIkMqa2VRZvPz8/LCwMEKIpaUl1ziKlZ6ezjW7at26NdfhNS0t7dmzZ4QQIyMjrl0TtxFuqqYrV664uLhwidzJsrCwMDU1lX0MSfWuTQD5knxMKK4oaehSyBkDj1QFoV5DlwJqqB7u5AjcPx9fYOBepbKysjVr1syePZsb3hugjkgN3D8IAnf4/CBwrwEE7nKtHu7kaCoDnzNlZWVe4weAulDZTE8AAAC1CKPKAAB8LDoUPVfvDgAAUBcQuAMA1Jy3t7eLiws3fry6urrUUf8BAABqBQJ3AICae/78ub+/f2FhoaGhoY+PD2rcAQCg7qCNOwBAzXl6enbs2NHIyKhTp050UCAAAIC6gMAdAKDmrK2t2Rm+AAAA6g6aygAAAAAAyAEE7gAAAAAAcgATMH0+amsCJgCof7g24fODCZhqABMwybV6uJOjxh0AAAAAQA4gcAcAAAAAkAMI3AEAAAAA5AACdwAAAAAAOYDAHQAAAABADmACJuDLj1id93hF/e9Xw25Z446r63+/AAAAAHIBNe4AAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAe+rMD9+vXr1tbWhw8flpFHJBKdOXNmzJgx9vb2Dg4O33333ZUrV2Tkj4uLW7p0qbOzc4cOHQYMGLBu3brk5GQZ+YODg6dOndqrV69OnToNGzbs0KFDJSWYWA4AAAAAqvBldU49evRodHR0RkZGZRlycnLc3d1v3bpFU0JCQo4dOzZy5Ehvb28VFRVe/pMnT06ePLmgoICmXL9+fevWrSdPnuzfv7/k9ufOnbt161b65+PHjy9cuLB9+3Y/Pz8zM7OafzEAAAAA+Nx9QYH7/fv3T58+LTvPt99+y0Xtjo6OgwYNKiwsvHjx4pMnT86ePaulpXXgwAE2c2Bg4HfffVdWVqajozN8+PCOHTsGBQVduHAhIyPD3d394cOH1tbWbP5ff/2Vi9qtrKyGDRtmaGh4/fr1K1euPH36dMiQIQ8fPhQKhbX9pQEAAADgMyEQi8UNXYa6VVJSEhUVdfLkyf379+fm5hJCNm3aNG/ePMmcly5dGjJkGTxBXgAAIABJREFUCCFk/PjxR44c4RLLy8sHDhx448YNQsiTJ0/s7Oxofnt7+7CwMA0NjaCgoHbt2nGJN27cGDhwYHl5+fDhw8+dO0czp6amtmjRorCw0NbWNigoqHHjxlz6ypUrV61aRQjZvn37zJkzP+abCgQCQojsE1qdPBgOEqD+VefaBJAvyceE4gq0Bf0wBh6pCkK9hi4F1FA93Mk/8zbuR48eVVNT69Kly+bNm7moXYadO3cSQszNzffv308TlZSUTp8+zdWF79q1i6YHBweHhYURQrZu3UqjdkJI//79586dSwj5888/379/T9MPHTpUWFgoEAjOnj1Lo3ZCyMqVKx0cHOjeAQAAAACk+sybyigpKTVv3pz++fbt28pyFhQUBAQEEEI8PT15bdm1tbX79et3+fLly5cv08RLly4RQho1ajR27FjepkaMGLFx40axWHzlypUpU6aw+Xv27Nm6dWtefnd39wcPHrx69SomJqZVq1Yf/i2/aOHh4bdu3UpISMjKytLT0zMzMxs4cGDLli0/fsuJiYlVPuyxjI2N2UcyVnFxcVxcnNRVGhoaxsbG3DM6AAAAgAyfeeA+duxYNrCWER5FRESUlpYSQnr16iW5dvDgwZcvX05MTExISDA2NiaEhIaGEkI6duyopqbGy9ylSxcDA4OUlJSHDx9ygbtIJHr06JGMjS9YsIAQ8vDhQwTu1efj47N06dI3b95IrrK3t9+0aVOfPn0+ZvszZ8709fX9oPJ4eHhERUUdP36cENKjR4+hQ4dyq6Kiouzt7Sv7oJqaWqtWraZNm/b9998rKHzmL8EAAACgxhAl/M/Lly+5BSsrK8m1bdu25RZiYmLY/FIzCwQCGxsbNvO7d++Ki4sry9+6dWslJSU2P8hWVFTk5ub27bffclG7vr6+s7PzyJEjHR0dNTQ0CCGhoaF9+/adN29e/bcYfvHixcaNGzdu3Hjt2rVqfqSwsDAiImL69On29vZ5eXl1WjwAAACQX595jXv10cHXmzVrJrnWyMiIW0hKSmLzS81M8/MyV5ZfQUHBwMAgISGB5q9MeHi47AyEEO69wWesrKzMxcXl7t27hJC2bdtu3rz5q6++UlRU5NYWFxefOnVq8eLFSUlJW7Zsyc3NZXssfJBdu3Zt2LCBl2hnZ1dYWDho0KDffvuNt8rQ0LDKba5bt27UqFFsSmJiYnh4+MaNG1NSUh49erRw4cI9e/bUrMDwGfjsr18AkK20tFRBAfcBqBQC9//hxmJXUFCQOiajqqoqm62iooKbNUmynQybn47vTheqmb8yXbp0qeJrMA8Jn6vly5dzUfuwYcNOnTrF65AgFArHjx8/aNCgAQMGPH78+MCBA7179/b09KzBjqQG4lxrKw0Njcqa0dvb23MzfLVp00ZyrZ6eHu+DLVu27N2797hx49q3b5+UlLR3795FixaZmprWoMDwGfjsr1/4oihjlKQPl5qaSlQqGroU8OlC4P4/XEUXrbjl4ZqyEEK4Fi+0VoymS83PZa5B/sp07txZxlquPl5ylqjPSWxs7ObNmwkh5ubmJ0+erOzL6urqXrx4sVWrVsXFxXPnzh05cmS9jZFvZmbm5eX1oZ/S1dVdsGDB7NmzxWLx06dPEbh/sT7v6xe+NAjba0BFRYXgPgCVQ+D+P1xdeFlZmVgsluzDSkNqrmpcVVVVIBCIxWKu3l0Sl5/W09OK9mrmrww3AGVluGJXp8GG/Nq2bVtFRQUhZM2aNY0aNZKR08TEZNasWRs3bkxNTT169OjUqVO59AMHDpSUlPTr14/rh8C6e/duZGRk8+bN3dzcalzCsrKyffv2EUKGDh36QfF3+/btuYXnz58PGjSoxgUAufZ5X7/wpUkWCBC7fyh9fX2M4w4yoHPq/9CB/KSOAEgTaTauE2RlwwVy6TSz7I1L5ofKcHNatW/f/ttvv60y8+LFi7W0tAgh7OAw8+fPnzFjRlBQkGT+06dPz5gx4yMH1C8pKZkxY8aMGTP+/vvvD/pgSkoKt9CkSZOPKQAAAAB8rhC4/0+LFi24BXbWJCohIYGXjVuQmpnm52WuLH9RUVFmZiabDaSKjo5OTEwkhLi7u1dn4PMmTZo4OzsTQgIDAz/9Pn9nz57lFtjZeQEAAAAoBO7/Q7sSRkZGSq7lEgUCgbW1NZtfauaKiopnz56x29TW1jYwMKgs/9OnT3llAKmio6O5BRljovNwOQsLC+mj1yeosLBw8eLF3MsEKyurDh06NHSJAAAA4FOENu7/Y21tbWxsnJCQcPPmzW+++Ya39ubNm4SQrl270tYs/fr1O3XqVExMTHx8PDs5KyHk4cOH+fn5hJD+/fvTxH79+vn4+Ny6dUty19zGlZWVe/fuXavf6XOTnp7OLejr61fzIzRnWlrap/BCY/fu3dwculRycvLz58+5tlICgeDAgQPonggAAABSocb9/40YMYIQcvr06dTUVDY9Ojr6xo0bhJDRo0fTxKFDh6qoqIjF4t27d/O2s2vXLkKIrq6uk5MTTeRG73716tXVq1fZzCUlJQcPHiSEuLi4aGpq1u43+szQrr2yu6WyaM5PpKnMkydPLvxbSEgIF7Wbmpr6+/vj4Q0AAAAqgxr3/zd//vx9+/bl5eWNGzfu7NmzXOV6SkrKN998IxKJ9PX16cgkhBBDQ8MpU6bs3r1727ZtvXr1Gjx4MJd+4MCBEydOEEIWLFjAxpfDhg1r167d06dPp06deuPGjVatWhFCSktLp02bFhsbKxAIli1bVq/fVg7p6OhwC1yXgOrIyMjgFnR1deukTB+oV69erVu3ZlMEAoGZmVmHDh369u3L9XgGAAAAkAqB+/8zMTHZs2fPpEmTrl27pq+v37t37/z8/IcPH5aXlyspKZ04cUJdXZ3Nv379+oCAgGfPng0ZMsTCwsLOzi44OJibP6V///6zZ89mMwsEAm9v7x49erx7965169adO3fW09O7d+8e16jm559/7tq1a31+WXlE56998eJFnz59qvORFy9eEEIEAkFlc9zWs/Hjx0+ePLmhSwEAAAByCU1l/mXixImnT59u1qxZcXHxtWvXgoKCysvLraysbt68yY1PwmrcuPGdO3dGjRqloKAQGxt7/vz55ORkFRWVH3/88eLFi5JzLbVr1y4wMLBTp05isTgsLOzKlSv5+fna2tp79uz5+eef6+sryrEuXbpw8yhxvQKqg8vZoUOH6oyxKMYkfwAAAPAJ+7Jq3KsTmY0aNcrd3T0wMDA2NlZRUdHa2lrGGCY6OjqnT59+//59SEhIZmamgYFBr169mjZtWln+9u3bh4eHR0ZGPnv2rKioyNTU1NHRsfottr9wQqGwT58+/v7+vr6+T58+bdeunWSeiIiI0tJS7pT5+/tzA9HQhkyyxcfH126BAQAAAGrRlxW4V5OiomLv3r2r303QxMTExMSk+ttv3749nSYTPsjs2bP9/f1FItGcOXOuX7/OW5uenu7m5paSknLw4EF3d/f58+cTQlRVVX/66SdeTqmPcHS4SQAAAIBPEJrKgDz5+uuv+/btSwi5cePGli1beGs1NDR69epVXFzs6enZtWtXboD8hQsXssNHKioqEkLi4uJ4n71+/frLly/rsuwAAAAAHwWBO8iZEydOcLNZzZs3z9PT8927d3SVUCg8duwYN/JmVFQUIcTJyYk3XI+VlRUh5PDhw3TAGULImzdvZsyYUT/lBwAAAKgZNJUBOdOsWbN79+65urq+evXK29vbx8ena9euNjY2TZo0SU5ODgkJiY2NpZnfv38fExND57slhIwYMeLBgwdJSUm2trZTpkwxMjKKiIg4depUYWHhtGnT9u7d2xDfCQAAAKBqCNxB/rRq1So0NHTdunU7d+4sLi5+8ODBgwcP2AyWlpYTJ048dOhQTExMt27dfHx8XF1duVVz5sx58uSJt7d3cnLy6tWruUQdHZ3z589nZGQgcAcAAIBPlgBD4H02BAIBqWrknOrkyY9Ynfd4Re2WrTo07JY17rj6gz6SnZ199erVu3fvJiYm5ufnGxkZmZubOzg4uLi4KCgoZGRkjBw5MiAgYNCgQZcuXWI/GBwcfP369devXwuFQjs7Ow8Pj6ZNm6akpDx//lxbW9vOzq6yPd67d6+iokJfX9/GxkZqhoqKinv37hFC2rdvT8cXysvLCw8PJ4S0bt36ExlRHj411bk2AeRL8jGhuKKkoUshZww8UhWEeg1dCqiheriTI3D/fHyBgXuVysrK1qxZM3v2bC0trdrdMkDtQuAOnx8E7jWAwF2u1cOdHE1l4HOmrKy8atWqhi4FAAAAQC3AqDIAAAAAAHIAgTsAAAAAgBxA4A4AAAAAIAcQuAMAAAAAyAEE7gAAAAAAcgCBOwAAAACAHEDgDgAAAAAgBzAB0+ejtiZgAoD6h2sTPj+YgKkGMAGTXKuHOzlq3AEAAAAA5AACdwAAAAAAOYDAHQAAAABADiBwBwAAAACQAwjcAQAAAADkgFJDFwA+OaufXF/x2L/+97vMrv/qTi71v18AAAAAuYAadwAAAAAAOYDAHQAAAABADiBwBwAAAACQAwjcAQAAAADkAAJ3AAAAAAA5gMAdAAAAAEAOIHAHAAAAAJADCNwBAAAAAOQAJmCCz8fLly9FIpGRkZGmpmZDlwUAAACgliFwh89H+/btS0pKjh075unp2dBlqX2//fZbcnIyIWTFihVqamoNXRwAAACobwjcAeTDwYMHnz17RgiZN28eAncAAIAvENq4AwAAAADIAdS4A8iHtWvXZmVlEUIaN27c0GUBAACABoDAHUA+uLm5NXQRAAAAoCEhcAf5k5GRsW/fvrt376alpenp6fXt23fSpEl6enpsnufPn9+8eZMQ4uXlpaGhwdtCSkrKmTNnCCFDhgwxMzMrKyvbt28fIWTChAnq6uo3b978448/Yv6PvfsMiOro+z4+C1IENCKiiA1RsYEFe0OxXdgriC32aDTGRI0mMbFEEjWWK8Uk9hJj72A02BUVRRSDqBgFQQRBighIh31enPveZ+8FEVH3cMz382qZM3vOfxddfgxzZu7dMzY2rl+//vjx49u0aVOwDLVa7e3t7ePjExUVlZiYWLly5Zo1a7q7u3fr1k2n54EDB2JiYnr27Ong4ODr67t169awsDBjY+MWLVrMmDGjdu3aQojw8PBff/316tWrmZmZ9vb2gwcPdnd31z7JxYsXg4KCbGxshg4dqlPGzp07Dx8+HBERYWpq2rRp0/Hjxzdr1uz69euXLl1q3Lixq6urdudDhw75+vqGhYWlpaXZ2dk5OztPmTKFUXwAAEo/lVqtlrsGvBkqlUoIUfQ3tDh9Ft84MT/I983WVhxfNe2+2Nntpd38/f2HDRsWFRWl3WhtbX3gwIHu3btrVpW5d++eg4ODEGLHjh3Dhw/XOcnSpUu/+OILExOT2NjYChUqpKWlScn14cOHS5cu/fXXX3X6f/3119988412S0JCgqura0hISMEKXV1djx07ZmJiomnp2LHjxYsXf//9d39//99++027c4UKFc6cORMeHj569Oj09HTtQx9++KF2JbNnz165cmWbNm0uX76saUxOTh49evSRI0e0n2hgYPDll1+am5t/8cUXEyZM2LBhg9T++PHjPn36BAUF6RRcqVKlQ4cOdejQoeBrgd4U5/8moCyx20zVeVlyV6EwVTyfGJhav7wfSiU9fJIz4g4lCQwM7Ny5c05OjoGBweDBgzt06JCWlvbXX39dvHhxwIABubm5mp716tVr1qzZjRs39u/fXzC4b9++XQjRv3//ChUqaLcvWbLkt99+6927t7u7e61ata5du7Z06dLExEQvL69evXq1a9dO03PYsGEhISEGBgYDBgzo2rWrpaXl48eP9+7dGxAQcObMma+++mr58uU6F12xYkVwcPCwYcM8PDxMTU137Nixffv25OTk/v37x8bGVqxY8dtvv23atGloaOiiRYvi4uLWrFkzefLkpk2bvujdyMnJ6dChw+3bt4UQ7dq169mzp5WV1YULF/bt2+fl5WVvb6/dOS8vb+jQoUFBQUZGRp6enh06dDA2Nr527dqmTZsSEhJGjBgREhLCuDsAAKUZwR1K8vnnn+fk5Jibm+/atatv375S47x587744otly5bpdPbw8Lhx48axY8fS09O1108MDg6WRspHjx6t85Tffvttzpw5mlO5urq6urq2adMmLy/v5MmTmuD+4MGD06dPCyGWLFkyZ84czdNnzZrVr1+/P//888SJEwWLDw4OXrRo0fz586Uve/funZKSIs20sbW1vXbtmo2NjXTR5s2bt2vXTq1WBwQEFBHcN27cKKX2efPmLV68WPpFf/r06YcOHRo+fHh4eLjO1S9duiSE2LBhw/vvvy81jhs3rnPnzh4eHg8fPvT39+/Zs+eLrgUAAGTHcpBQDD8/P2na+qxZszSpXQihUqmWLl3q7Oys09/Dw0MIkZ6efvToUe12abjd2tq6V69eOk+pVavW4sWLtVtatGjh6OgohNDOwcHBwVZWVpUrV/7www+1O6tUqkGDBgkhwsLCCtZvZ2f35ZdfardInYUQ3377rZTaJW3btrW1tRVC6MwI0padne3l5SWE6Nixo5eXl5TaJQMHDpw6dapOf80MmYEDB2q3Dxw40N3dfciQIdp/rwAAAKUQwR2KIc3tNjEx+eSTTwoenTlzpk5LnTp1WrZsKYTYt2+fplG6lVMIMXz48DJldP/iNHLkSGNjY53G6tWri/87ZW3AgAEJCQlxcXEF55ZIKzYWOr+tffv2OlesWrWq9KBz5846naVD+fn5Bc8juXv3bnR0tBDi888/L3j0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Laxr167Xrl07cOBARESEzgR6b2/vK1euGBkZOTg4GBgY9O3bNzExccKECRs2bNDuVrZsWScnp9DQ0JKtLAQAQOkh21SZyMjIqVOnSqnd0tKyWbNmLV6smLs54t9gx44d0iLis2fPHjVqlPbmu6amptu2bXN3dxdCSKnd1dVVZ38faZHvzZs3a9ZFEUI8ePBg+vTpOhcyMzOT7onUmVKfnZ09a9asM2fOCCGysrKkcP86dFZUlCrMy8v773//q92+ZMkSadz6Lfnqq69cXV2FEIcOHerYsWPBa8XExEyYMGHy5MlCCEtLy927d0u/Y+gYMWKE9q9AQUFBUr5/7733pL+BaNN57QV9/vnnFSpUyM3N9fDwkHaJkoSHh8+bN08I0a1bN2l9d2dnZyHE3r17b9++rX2GmzdvHj9+XAjRrl07TePmzZs9PT09PT39/f2LLgAAgNJDthH3Y8eOZWRkGBgYbNy4ccyYMcXZwR4QQlStWtXPz693797379/fvn37zp07W7du3ahRo/feey82NvbKlSvaaw4+evTo3r17DRo00LQMGTLk8uXLjx8/dnR0nDRpkq2t7d9//7179+709PQpU6bozPD29PRcvnz5vn373NzcevXqZWpqeu/evT179kRFRXXq1MnPzy83N7d3797jx48fN25cCV6LtEmwn5+fm5tbpUqVfvzxRysrq3bt2lWtWvXx48dffvnl1atXXV1dnz17dvLkyXPnzvXr1y8kJOTBgwclffOKYmho6O3t7enp+eeffwYGBnbs2LFu3bpt27a1trZOSUm5c+fO5cuXpUV1atasefToUXt7+4InqVWr1r1791q0aNG5c2cHB4fg4OCLFy9K82oWL16svSlyoa+90Lfo22+/nTZt2tWrV52cnLp06eLg4BAZGXngwIHMzMz33ntP8y2bO3fuqVOnUlJSWrVq1b9/f3t7+9zc3Dt37hw7diw3N9fGxmbatGma0167dk1admbgwIHagR4AgNJMtuAupSsPDw/NjktAMdWrV+/q1avffffdzz//nJmZefny5cuXL2t3qFOnzvjx4zdu3Hjv3r02bdrs3Lmzd+/e0qGZM2feuHFj+/btsbGxixcvlhqtrKwOHjyYmJioE9y9vLxCQ0N9fHx8fX19fX2lxvLly2/ZsmXMmDHDhg3bs2fPhQsXzM3NSxbcx4wZs27d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"/>
```

---

## System Information

````julia
using InteractiveUtils

versioninfo()
````

````
Julia Version 1.11.1
Commit 8f5b7ca12ad (2024-10-16 10:53 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 32 × 13th Gen Intel(R) Core(TM) i9-13900KF
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, alderlake)
Threads: 16 default, 0 interactive, 8 GC (on 32 virtual cores)
Environment:
  JULIA_CONDAPKG_BACKEND = Null
  JULIA_PYTHONCALL_EXE = /home/alberto/.julia/dev/QuantumToolbox/benchmarks/package_comparison/pyenv/bin/python
  LD_LIBRARY_PATH = /usr/local/lib:
  JULIA_NUM_THREADS = 16

````

---

````julia
QuantumToolbox.about()
````

````

 QuantumToolbox.jl: Quantum Toolbox in Julia
≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
Copyright © QuTiP team 2022 and later.
Current admin team:
    Alberto Mercurio and Yi-Te Huang

Package information:
====================================
Julia              Ver. 1.11.1
QuantumToolbox     Ver. 0.21.5
SciMLOperators     Ver. 0.3.12
LinearSolve        Ver. 2.37.0
OrdinaryDiffEqCore Ver. 1.11.0

System information:
====================================
OS       : Linux (x86_64-linux-gnu)
CPU      : 32 × 13th Gen Intel(R) Core(TM) i9-13900KF
Memory   : 62.514 GB
WORD_SIZE: 64
LIBM     : libopenlibm
LLVM     : libLLVM-16.0.6 (ORCJIT, alderlake)
BLAS     : libopenblas64_.so (ilp64)
Threads  : 16 (on 32 virtual cores)


````

---

````julia
qutip.about()
````

````
Python: None
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

