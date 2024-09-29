# [Time Evolution and Quantum System Dynamics](@id doc:Time-Evolution-and-Quantum-System-Dynamics)

**Table of contents**

```@contents
Pages = [
    "intro.md",
    "solution.md",
    "sesolve.md",
    "mesolve.md",
    "mcsolve.md",
    "stochastic.md",
    "time_dependent.md",
]
Depth = 2:3
```

## [Introduction](@id doc-TE:Introduction)

Although in some cases, we want to find the stationary states of a quantum system, often we are interested in the dynamics: how the state of a system or an ensemble of systems evolves with time. `QuantumToolbox` provides many ways to model dynamics.

There are two kinds of quantum systems: open systems that interact with a larger environment and closed systems that do not. In a closed system, the state can be described by a state vector. When we are modeling an open system, or an ensemble of systems, the use of the density matrix is mandatory.

The following table lists the solvers provided by `QuantumToolbox` for dynamic quantum systems and the corresponding type of solution returned by the solver:

| **Equation** | **Function Call** | **Problem** | **Returned Solution** |
|:-------------|:------------------|:------------|:----------------------|
| Unitary evolution, Schrödinger equation | [`sesolve`](@ref) | [`sesolveProblem`](@ref) | [`TimeEvolutionSol`](@ref) |
| Lindblad master eqn. or Von Neuman eqn.| [`mesolve`](@ref) | [`mesolveProblem`](@ref) | [`TimeEvolutionSol`](@ref) |
| Monte Carlo evolution | [`mcsolve`](@ref) | [`mcsolveProblem`](@ref) [`mcsolveEnsembleProblem`](@ref) | [`TimeEvolutionMCSol`](@ref) |
| Stochastic Schrödinger equation | [`ssesolve`](@ref) | [`ssesolveProblem`](@ref) [`ssesolveEnsembleProblem`](@ref) | [`TimeEvolutionSSESol`](@ref) |

!!! note "Solving dynamics with pre-defined problems"
    `QuantumToolbox` provides two different methods to solve the dynamics. One can use the function calls listed above by either taking all the operators (like Hamiltonian and collapse operators, etc.) as inputs directly, or generating the `prob`lems by yourself and take it as an input of the function call, e.g., `sesolve(prob)`.
