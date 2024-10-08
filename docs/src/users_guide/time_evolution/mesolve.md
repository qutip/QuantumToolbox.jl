# [Lindblad Master Equation Solver](@id doc-TE:Lindblad-Master-Equation-Solver)

```@setup mesolve
using QuantumToolbox
```

## Von Neumann equation

While the evolution of the state vector in a closed quantum system is deterministic (as we discussed in the previous section: [Schrödinger Equation Solver](@ref doc-TE:Schrödinger-Equation-Solver)), open quantum systems are stochastic in nature. The effect of an environment on the system of interest is to induce stochastic transitions between energy levels, and to introduce uncertainty in the phase difference between states of the system. The state of an open quantum system is therefore described in terms of ensemble averaged states using the density matrix formalism. A density matrix ``\hat{\rho}`` describes a probability distribution of quantum states ``|\psi_n\rangle`` in a matrix representation, namely

```math
\hat{\rho} = \sum_n p_n |\psi_n\rangle\langle\psi_n|,
```

where ``p_n`` is the classical probability that the system is in the quantum state ``|\psi_n\rangle``. The time evolution of a density matrix ``\hat{\rho}`` is the topic of the remaining portions of this section.

The time evolution of the density matrix ``\hat{\rho}(t)`` under closed system dynamics is governed by the von Neumann equation: 

```math
\frac{d\hat{\rho}(t)}{dt} = -\frac{i}{\hbar}\left[\hat{H}, \hat{\rho}(t)\right],
```

where ``[\cdot, \cdot]`` represents the commutator. In `QuantumToolbox`, given a Hamiltonian, we can calculate the unitary (non-dissipative) time-evolution of an arbitrary initial state using the `QuantumToolbox` time evolution problem [`mesolveProblem`](@ref) or directly call the function [`mesolve`](@ref). It evolves the density matrix ``\hat{\rho}(t)`` and evaluates the expectation values for a set of operators `e_ops` at each given time points, using an ordinary differential equation solver provided by the powerful julia package [`DifferentialEquation.jl`](https://docs.sciml.ai/DiffEqDocs/stable/).

```@example mesolve
H = 0.5 * sigmax()
state0 = basis(2, 0)         # state vector
state0 = ket2dm(basis(2, 0)) # density matrix
tlist = LinRange(0.0, 10.0, 20)

sol = mesolve(H, state0, tlist, e_ops = [sigmaz()])
```

!!! note "Type of initial state"
    The initial state `state0` here can be given as a state vector ``|\psi(0)\rangle`` (in the type of [`Ket`](@ref)) or a density matrix ``\hat{\rho}(0)`` (in the type of [`Operator`](@ref)). If it is given as a [`Ket`](@ref), it will be transformed to density matrix ``\hat{\rho}(0) = |\psi(0)\rangle\langle\psi(0)|`` internally in [`mesolve`](@ref).

The function returns [`TimeEvolutionSol`](@ref), as described in the previous section [Time Evolution Solutions](@ref doc-TE:Time-Evolution-Solutions). The stored `states` will always be in the type of [`Operator`](@ref) (density matrix).

```@example mesolve
sol.states
```

Here, only the final state is stored because the `states` will be saved depend on the keyword argument `saveat` in `kwargs`. If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). 

One can also specify `e_ops` and `saveat` separately:

```@example mesolve
tlist = [0, 5, 10]
sol = mesolve(H, state0, tlist, e_ops = [sigmay()], saveat = tlist)
```

```@example mesolve
sol.expect
```

```@example mesolve
sol.states
```

## The Lindblad master equation