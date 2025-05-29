# [Lindblad Master Equation Solver](@id doc-TE:Lindblad-Master-Equation-Solver)

```@setup mesolve
using QuantumToolbox

using CairoMakie
CairoMakie.enable_only_mime!(MIME"image/svg+xml"())
```

## [Von Neumann equation](@id doc-TE:Von-Neumann-equation)

While the evolution of the state vector in a closed quantum system is deterministic (as we discussed in the previous section: [Schrödinger Equation Solver](@ref doc-TE:Schrödinger-Equation-Solver)), open quantum systems are stochastic in nature. The effect of an environment on the system of interest is to induce stochastic transitions between energy levels, and to introduce uncertainty in the phase difference between states of the system. The state of an open quantum system is therefore described in terms of ensemble averaged states using the density matrix formalism. A density matrix ``\hat{\rho}`` describes a probability distribution of quantum states ``|\psi_n\rangle`` in a matrix representation, namely

```math
\hat{\rho} = \sum_n p_n |\psi_n\rangle\langle\psi_n|,
```

where ``p_n`` is the classical probability that the system is in the quantum state ``|\psi_n\rangle``. The time evolution of a density matrix ``\hat{\rho}`` is the topic of the remaining portions of this section.

The time evolution of the density matrix ``\hat{\rho}(t)`` under closed system dynamics is governed by the von Neumann equation: 

```math
\begin{equation}
\frac{d}{dt}\hat{\rho}(t) = -\frac{i}{\hbar}\left[\hat{H}, \hat{\rho}(t)\right],
\end{equation}
```

where ``[\cdot, \cdot]`` represents the commutator. The above equation is equivalent to the Schrödinger equation described in the [previous section](@ref doc-TE:Schrödinger-Equation-Solver) under the density matrix formalism.

In `QuantumToolbox`, given a Hamiltonian, we can calculate the unitary (non-dissipative) time-evolution of an arbitrary initial state using the `QuantumToolbox` time evolution problem [`mesolveProblem`](@ref) or directly call the function [`mesolve`](@ref). It evolves the density matrix ``\hat{\rho}(t)`` and evaluates the expectation values for a set of operators `e_ops` at each given time points, using an ordinary differential equation solver provided by the powerful julia package [`DifferentialEquation.jl`](https://docs.sciml.ai/DiffEqDocs/stable/).

```@example mesolve
H = 0.5 * sigmax()
state0 = basis(2, 0) # state vector
tlist = LinRange(0.0, 10.0, 20)

sol = mesolve(H, state0, tlist, e_ops = [sigmaz()])
```

!!! note "Use sesolve for improved efficiency"
    Here, if the Hamiltonian `H` is given as an [`Operator`](@ref), and the initial state `state0` is given as a state vector ``|\psi(0)\rangle`` (in the type of [`Ket`](@ref)), it will automatically call [`sesolve`](@ref) for improved efficiency.

The function returns [`TimeEvolutionSol`](@ref), as described in the previous section [Time Evolution Solutions](@ref doc-TE:Time-Evolution-Solutions).

```@example mesolve
sol.states
```

Here, only the final state is stored because the `states` will be saved depend on the keyword argument `saveat` in `kwargs`. If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). 

One can also specify `e_ops` and `saveat` separately:

```@example mesolve
tlist = [0, 5, 10]
state0 = ket2dm(basis(2, 0))  # density matrix
sol = mesolve(H, state0, tlist, e_ops = [sigmay()], saveat = tlist)
```

```@example mesolve
sol.expect
```

```@example mesolve
sol.states
```

Note that when the initial state `state0` is given as a density matrix ``|\psi(0)\rangle\langle\psi(0)|`` (in the type of [`Operator`](@ref)), the stored `states` will also be in the type of [`Operator`](@ref) (density matrix).

## [The Lindblad master equation](@id doc-TE:The-Lindblad-master-equation)

The standard approach for deriving the equations of motion for a system interacting with its environment is to expand the scope of the system to include the environment. The combined quantum system is then closed, and its evolution is also governed by the von Neumann equation

```math
\begin{equation}
\frac{d}{dt}\hat{\rho}_{\textrm{tot}}(t) = -\frac{i}{\hbar}\left[\hat{H}_{\textrm{tot}}, \hat{\rho}_{\textrm{tot}}(t)\right].
\end{equation}
```

Here, the total Hamiltonian

```math
\hat{H}_{\textrm{tot}} = \hat{H}_{\textrm{sys}} + \hat{H}_{\textrm{env}} + \hat{H}_{\textrm{int}},
```

includes the original system Hamiltonian ``\hat{H}_{\textrm{sys}}``, the Hamiltonian for the environment ``\hat{H}_{\textrm{env}}``, and a term representing the interaction between the system and its environment ``\hat{H}_{\textrm{int}}``. Since we are only interested in the dynamics of the system, we can, perform a partial trace over the environmental degrees of freedom, and thereby obtain a master equation for the motion of the original system density matrix ``\hat{\rho}_{\textrm{sys}}(t)=\textrm{Tr}_{\textrm{env}}[\hat{\rho}_{\textrm{tot}}(t)]``. The most general trace-preserving and completely positive form of this evolution is the Lindblad master equation for the reduced density matrix, namely

```math
\begin{equation}
\frac{d}{dt}\hat{\rho}_{\textrm{sys}}(t) = -\frac{i}{\hbar}\left[\hat{H}_{\textrm{sys}}, \hat{\rho}_{\textrm{sys}}(t)\right] + \sum_n \hat{C}_n \hat{\rho}_{\textrm{sys}}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}_{\textrm{sys}}(t) - \frac{1}{2} \hat{\rho}_{\textrm{sys}}(t) \hat{C}_n^\dagger \hat{C}_n
\end{equation}
```

where ``\hat{C}_n \equiv \sqrt{\gamma_n}\hat{A}_n`` are the collapse operators, ``\hat{A}_n`` are the operators acting on the system in ``\hat{H}_{\textrm{int}}`` which describes the system-environment interaction, and ``\gamma_n`` are the corresponding rates. The derivation of Lindblad master equation may be found in several sources, and will not be reproduced here. Instead, we emphasize the approximations that are required to arrive at the above Lindblad master equation from physical arguments, and hence perform a calculation in `QuantumToolbox`:

- **Separability:** At ``t = 0``, there are no correlations between the system and environment, such that the total density matrix can be written as a tensor product, namely ``\hat{\rho}_{\textrm{tot}}(0)=\hat{\rho}_{\textrm{sys}}(0)\otimes\hat{\rho}_{\textrm{env}}(0)``.
- **Born approximation:** Requires: (i) the state of the environment does not significantly change as a result of the interaction with the system; (ii) the system and the environment remain separable throughout the evolution. These assumptions are justified if the interaction is weak, and if the environment is much larger than the system. In summary, ``\hat{\rho}_{\textrm{tot}}(t)\approx\hat{\rho}_{\textrm{sys}}(t)\otimes\hat{\rho}_{\textrm{env}}(0)``.
- **Markov approximation:** The time-scale of decay for the environment ``\tau_{\textrm{env}}`` is much shorter than the smallest time-scale of the system dynamics, i.e., ``\tau_{\textrm{sys}}\gg\tau_{\textrm{env}}``. This approximation is often deemed a “short-memory environment” as it requires the environmental correlation functions decay in a fast time-scale compared to those of the system.
- **Secular approximation:** Stipulates that elements in the master equation corresponding to transition frequencies satisfy ``|\omega_{ab}-\omega_{cd}| \ll 1/\tau_{\textrm{sys}}``, i.e., all fast rotating terms in the interaction picture can be neglected. It also ignores terms that lead to a small renormalization of the system energy levels. This approximation is not strictly necessary for all master-equation formalisms (e.g., the Block-Redfield master equation), but it is required for arriving at the Lindblad form in the above equation which is used in [`mesolve`](@ref).

For systems with environments satisfying the conditions outlined above, the Lindblad master equation governs the time-evolution of the system density matrix, giving an ensemble average of the system dynamics. In order to ensure that these approximations are not violated, it is important that the decay rates ``\gamma_n`` be smaller than the minimum energy splitting in the system Hamiltonian. Situations that demand special attention therefore include, for example, systems strongly coupled to their environment, and systems with degenerate or nearly degenerate energy levels.

What is new in the master equation compared to the Schrödinger equation (or von Neumann equation) are processes that describe dissipation in the quantum system due to its interaction with an environment. For example, evolution that includes incoherent processes such as relaxation and dephasing. These environmental interactions are defined by the operators ``\hat{A}_n`` through which the system couples to the environment, and rates ``\gamma_n`` that describe the strength of the processes.

In `QuantumToolbox`, the function [`mesolve`](@ref) can also be used for solving the master equation. This is done by passing a list of collapse operators (`c_ops`) as the fourth argument of the [`mesolve`](@ref) function in order to define the dissipation processes of the Lindblad master equation. As we mentioned above, each collapse operator ``\hat{C}_n`` is the product of ``\sqrt{\gamma_n}`` (the square root of the rate) and ``\hat{A}_n`` (operator which describes the dissipation process). 

Furthermore, `QuantumToolbox` solves the master equation in the [`SuperOperator`](@ref) formalism. That is, a Liouvillian [`SuperOperator`](@ref) will be generated internally in [`mesolve`](@ref) by the input system Hamiltonian ``\hat{H}_{\textrm{sys}}`` and the collapse operators ``\hat{C}_n``. One can also generate the Liouvillian [`SuperOperator`](@ref) manually for special purposes, and pass it as the first argument of the [`mesolve`](@ref) function. To do so, it is useful to read the section [Superoperators and Vectorized Operators](@ref doc:Superoperators-and-Vectorized-Operators), and also the docstrings of the following functions:
- [`spre`](@ref)
- [`spost`](@ref)
- [`sprepost`](@ref)
- [`liouvillian`](@ref)
- [`lindblad_dissipator`](@ref)

## [Example: Dissipative Spin dynamics](@id doc-TE:Example:Dissipative-Spin-dynamics)

Using the example with the dynamics of spin-``\frac{1}{2}`` from the previous section ([Schrödinger Equation Solver](@ref doc-TE:Schrödinger-Equation-Solver)), we can easily add a relaxation process (describing the dissipation of energy from the spin to the environment), by adding `[sqrt(γ) * sigmax()]` in the fourth parameter of the [`mesolve`](@ref) function.

```@example mesolve
H = 2 * π * 0.1 * sigmax()
ψ0 = basis(2, 0) # spin-up
tlist = LinRange(0.0, 10.0, 100)

γ = 0.05
c_ops = [sqrt(γ) * sigmax()]

sol = mesolve(H, ψ0, tlist, c_ops, e_ops = [sigmaz(), sigmay()])
```

We can therefore plot the expectation values:

```@example mesolve
times = sol.times
expt_z = real(sol.expect[1,:])
expt_y = real(sol.expect[2,:])

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Expectation values")
lines!(ax, times, expt_z, label = L"\langle\hat{\sigma}_z\rangle", linestyle = :solid)
lines!(ax, times, expt_y, label = L"\langle\hat{\sigma}_y\rangle", linestyle = :dash)

axislegend(ax, position = :rt)

fig
```

## [Example: Harmonic oscillator in thermal bath](@id doc-TE:Example:Harmonic-oscillator-in-thermal-bath)

Consider a harmonic oscillator (single-mode cavity) couples to a thermal bath. If the single-mode cavity initially is in a `10`-photon [`fock`](@ref) state, the dynamics is calculated with the following code:

```@example mesolve
# Define parameters
N = 20  # number of basis states to consider
a = destroy(N)
H = a' * a
ψ0 = fock(N, 10)  # initial state
κ = 0.1  # coupling to oscillator
n_th = 2 # temperature with average of 2 excitations
tlist = LinRange(0, 50, 100)

# collapse operators 
c_ops = [
    sqrt(κ * (n_th + 1)) * a, # emission
    sqrt(κ *  n_th     ) * a' # absorption
]

# find expectation value for particle number
e_ops = [a' * a]

sol = mesolve(H, ψ0, tlist, c_ops, e_ops=e_ops)
Num = real(sol.expect[1, :])

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], 
    title = L"Decay of Fock state $|10\rangle$ in a thermal environment with $\langle n\rangle=2$",
    xlabel = "Time", 
    ylabel = "Number of excitations",
)
lines!(ax, tlist, Num)

fig
```

## [Example: Two-level atom coupled to dissipative single-mode cavity](@id doc-TE:Example:Two-level-atom-coupled-to-dissipative-single-mode-cavity)

Consider a two-level atom coupled to a dissipative single-mode cavity through a dipole-type interaction, which supports a coherent exchange of quanta between the two systems. If the atom initially is in its ground state and the cavity in a `5`-photon [`fock`](@ref) state, the dynamics is calculated with the following code:

!!! note "Generate Liouviilian superoperator manually"
    In this example, we demonstrate how to generate the Liouvillian [`SuperOperator`](@ref) manually and pass it as the first argument of the [`mesolve`](@ref) function.

```@example mesolve
# two-level atom
σm = tensor(destroy(2), qeye(10))
H_a = 2 * π * σm' * σm

# dissipative single-mode cavity
γ = 0.1 # dissipation rate
a  = tensor(qeye(2), destroy(10))
H_c = 2 * π * a' * a
c_ops = [sqrt(γ) * a]

# coupling between two-level atom and single-mode cavity
g = 0.25 # coupling strength
H_I = 2 * π * g * (σm * a' + σm' * a)

ψ0 = tensor(basis(2,0), fock(10, 5)) # initial state
tlist = LinRange(0.0, 10.0, 200)

# generate Liouvillian superoperator manually
L = liouvillian(H_a + H_c + H_I, c_ops)
sol = mesolve(L, ψ0, tlist, e_ops=[σm' * σm, a' * a])

times = sol.times

# expectation value of Number operator
N_atom   = real(sol.expect[1,:])
N_cavity = real(sol.expect[2,:])

# plot by CairoMakie.jl
fig = Figure(size = (500, 350))
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Expectation values")
lines!(ax, times, N_atom, label = "atom excitation probability", linestyle = :solid)
lines!(ax, times, N_cavity, label = "cavity photon number", linestyle = :dash)

axislegend(ax, position = :rt)

fig
```
