# [Monte Carlo Solver](@id doc-TE:Monte-Carlo-Solver)

## [Monte Carlo wave-function](@id doc-TE:Monte-Carlo-wave-function)

Where as the density matrix formalism describes the ensemble average over many identical realizations of a quantum system, the Monte Carlo (MC), or quantum-jump approach to wave function evolution, allows for simulating an individual realization of the system dynamics. Here, the environment is continuously monitored, resulting in a series of quantum jumps in the system wave function, conditioned on the increase in information gained about the state of the system via the environmental measurements. In general, this evolution is governed by the Schrödinger equation with a non-Hermitian effective Hamiltonian

```math
\hat{H}_{\textrm{eff}} = \hat{H} - \frac{i}{2} \sum_{n=1}^N \hat{C}_n^\dagger \hat{C}_n.
```

where ``\hat{H}`` is the system Hamiltonian and ``\hat{C}_n`` are collapse (jump) operators (assume ``N`` is the total number of collapse operators). Each collapse operator corresponds to a separate irreversible process with rate ``\gamma_n``. Here, the strictly negative non-Hermitian portion of the above equation gives rise to a reduction in the norm of the wave function, that to first-order in a small time ``\delta t``, is given by

```math
\langle \psi(t + \delta t) | \psi(t + \delta t) \rangle = 1 - \delta p,
```

where

```math
\delta p = \delta t \sum_{n=1}^N \langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle,
```

and ``\delta t`` is such that ``\delta p \ll 1``. With a probability of remaining in the state ``| \psi(t + \delta t) \rangle`` given by ``1 - \delta p``, the corresponding quantum jump probability is thus ``\delta p``. If the environmental measurements register a quantum jump, say via the emission of a photon into the environment, or a change in the spin of a quantum dot, the wave function undergoes a jump into a state defined by projecting ``| \psi(t) \rangle`` using the collapse operator ``\hat{C}_n`` corresponding to the measurement

```math
| \psi(t+\delta t) \rangle = \frac{\hat{C}_n |\psi(t)\rangle}{ \sqrt{\langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle} }.
```

If more than a single collapse operator is present in ``\hat{H}_{\textrm{eff}}``, the probability of collapse due to the ``n``-th operator ``\hat{C}_n`` is given by

```math
P_n(t) = \frac{1}{\delta p}\langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle.
```

Note that the probability of all collapses should be normalized to unity for all time ``t``, namely

```math
\sum_{n=1}^N P_n(t) = 1 ~~~\forall~~t.
```

Evaluating the MC evolution to first-order in time is quite tedious. Instead, `QuantumToolbox.jl` provides the function [`mcsolve`](@ref) which uses the following algorithm to simulate a single realization of a quantum system. Starting from a pure state ``| \psi(0) \rangle``:

1. Choose two random numbers (``r_1`` and ``r_2``) between 0 and 1, where ``r_1`` represents the probability that a quantum jump occurs and  ``r_2`` is used to select which collapse operator was responsible for the jump.
1. Integrate the Schrödinger equation with respect to the effective Hamiltonian ``\hat{H}_{\textrm{eff}}`` until a time ``\tau`` such that the norm of the wave function satisfies ``\langle \psi(\tau) | \psi(\tau) \rangle = r_1``, at which point a jump occurs
1. The resultant jump projects the system at time ``\tau`` into one of the renormalized states ``| \psi(\tau + \delta t) \rangle``. The corresponding collapse operator ``\hat{C}_n`` is chosen such that ``\tilde{n} \leq N`` is the smallest integer satisfying ``\sum_{n=1}^{\tilde{n}} P_n(\tau) \geq r_2``.
1. Using the renormalized state from previous step as the new initial condition at time ``\tau`` and repeat the above procedure until the final simulation time is reached.

## [Example: Two-level atom coupled to dissipative single-mode cavity (MC)](@id doc-TE:Example:Two-level-atom-coupled-to-dissipative-single-mode-cavity-(MC))

In `QuantumToolbox.jl`, Monte Carlo evolution is implemented with the [`mcsolve`](@ref) function. It takes nearly the same arguments as the [`mesolve`](@ref) function for [Lindblad master equation evolution](@ref doc-TE:Lindblad-Master-Equation-Solver), except that the initial state must be a [`Ket`](@ref) vector, as oppose to a density matrix, and there is an optional keyword argument `ntraj` that defines the number of stochastic trajectories to be simulated. By default, `ntraj=500` indicating that `500` Monte Carlo trajectories will be performed.

To illustrate the use of the Monte Carlo evolution of quantum systems in `QuantumToolbox.jl`, let’s again consider the case of a two-level atom coupled to a leaky cavity. The only differences to the master equation treatment is that in this case we invoke the [`mcsolve`](@ref) function instead of [`mesolve`](@ref)

```@setup mcsolve
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

```@example mcsolve
times = LinRange(0.0, 10.0, 200)

ψ0 = tensor(fock(2, 0), fock(10, 8))
a  = tensor(qeye(2), destroy(10))
σm = tensor(destroy(2), qeye(10))
H = 2 * π * a' * a + 2 * π * σm' * σm + 2 * π * 0.25 * (σm * a' + σm' * a)

c_ops = [sqrt(0.1) * a]
e_ops = [a' * a, σm' * σm]

sol_500 = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Expectation values",
    title = "Monte Carlo time evolution (500 trajectories)",
)
lines!(ax, times, real(sol_500.expect[1,:]), label = "cavity photon number", linestyle = :solid)
lines!(ax, times, real(sol_500.expect[2,:]), label = "atom excitation probability", linestyle = :dash)

axislegend(ax, position = :rt)

fig
```

The advantage of the Monte Carlo method over the master equation approach is that only the state vector ([`Ket`](@ref)) is required to be kept in the computers memory, as opposed to the entire density matrix. For large quantum system this becomes a significant advantage, and the Monte Carlo solver is therefore generally recommended for such systems. However, for small systems, the added overhead of averaging a large number of stochastic trajectories to obtain the open system dynamics, as well as starting the multiprocessing functionality, outweighs the benefit of the minor (in this case) memory saving. Master equation methods are therefore generally more efficient when Hilbert space sizes are on the order of a couple of hundred states or smaller.

We can also change the number of trajectories (`ntraj`). This can be used to explore the convergence of the Monte Carlo solver. For example, the following code plots the expectation values for `1`, `10` and `100` trajectories:

```@example mcsolve
e_ops = [a' * a]

sol_1   = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops, ntraj = 1)
sol_10  = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops, ntraj = 10)
sol_100 = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops, ntraj = 100)

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Expectation values",
    title = "Monte Carlo time evolution",
)
lines!(ax, times, real(sol_1.expect[1,:]), label = "1 trajectory", linestyle = :dashdot)
lines!(ax, times, real(sol_10.expect[1,:]), label = "10 trajectories", linestyle = :dash)
lines!(ax, times, real(sol_100.expect[1,:]), label = "100 trajectories", linestyle = :solid)

axislegend(ax, position = :rt)

fig
```

## [Running trajectories in parallel](@id doc-TE:Running-trajectories-in-parallel)

Monte Carlo evolutions often need hundreds of trajectories to obtain sufficient statistics. Since all trajectories are independent of each other, they can be computed in parallel. The keyword argument `ensemblealg` can specify how the multiple trajectories are handled. The common ensemble methods are:

- `EnsembleSerial()`: No parallelism
- `EnsembleThreads()`: **The default.** This uses multithreading.
- `EnsembleDistributed()`: This uses as many processors as you have Julia processes.
- `EnsembleSplitThreads()`: This uses multithreading on each process.

!!! note "Other Ensemble Algorithms"
    See the [documentation of `DifferentialEquations.jl`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/) for more details. Also, see Julia's documentation for more details about multithreading and adding more processes.

```julia
sol_serial   = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops, ensemblealg=EnsembleSerial())
sol_parallel = mcsolve(H, ψ0, times, c_ops, e_ops=e_ops, ensemblealg=EnsembleThreads());
```

!!! tip "Parallelization on a Cluster"
    See the section [Intensive parallelization on a Cluster](@ref doc:Intensive-parallelization-on-a-Cluster) for more details.
