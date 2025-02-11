# [Time Evolution and Quantum System Dynamics](@id doc:Time-Evolution-and-Quantum-System-Dynamics)

**Table of contents**

- [Introduction](@ref doc-TE:Introduction)
- [Time Evolution Solutions](@ref doc-TE:Time-Evolution-Solutions)
    - [Solution](@ref doc-TE:Solution)
    - [Accessing data in solutions](@ref doc-TE:Accessing-data-in-solutions)
    - [Multiple trajectories solution](@ref doc-TE:Multiple-trajectories-solution)
- [Schrödinger Equation Solver](@ref doc-TE:Schrödinger-Equation-Solver)
    - [Unitary evolution](@ref doc-TE:Unitary-evolution)
    - [Example: Spin dynamics](@ref doc-TE:Example:Spin-dynamics)
- [Lindblad Master Equation Solver](@ref doc-TE:Lindblad-Master-Equation-Solver)
    - [Von Neumann equation](@ref doc-TE:Von-Neumann-equation)
    - [The Lindblad master equation](@ref doc-TE:The-Lindblad-master-equation)
    - [Example: Dissipative Spin dynamics](@ref doc-TE:Example:Dissipative-Spin-dynamics)
    - [Example: Harmonic oscillator in thermal bath](@ref doc-TE:Example:Harmonic-oscillator-in-thermal-bath)
    - [Example: Two-level atom coupled to dissipative single-mode cavity](@ref doc-TE:Example:Two-level-atom-coupled-to-dissipative-single-mode-cavity)
- [Monte Carlo Solver](@ref doc-TE:Monte-Carlo-Solver)
- [Stochastic Solver](@ref doc-TE:Stochastic-Solver)
    - [Stochastic Schrödinger equation](@ref doc-TE:Stochastic-Schrödinger-equation)
    - [Stochastic master equation](@ref doc-TE:Stochastic-master-equation)
    - [Example: Homodyne detection](@ref doc-TE:Example:Homodyne-detection)
- [Solving Problems with Time-dependent Hamiltonians](@ref doc-TE:Solving-Problems-with-Time-dependent-Hamiltonians)
    - [Generate QobjEvo](@ref doc-TE:Generate-QobjEvo)
    - [QobjEvo fields (attributes)](@ref doc-TE:QobjEvo-fields-(attributes))
    - [Using parameters](@ref doc-TE:Using-parameters)

# [Introduction](@id doc-TE:Introduction)

Although in some cases, we want to find the stationary states of a quantum system, often we are interested in the dynamics: how the state of a system or an ensemble of systems evolves with time. `QuantumToolbox` provides many ways to model dynamics.

There are two kinds of quantum systems: open systems that interact with a larger environment and closed systems that do not. In a closed system, the state can be described by a state vector. When we are modeling an open system, or an ensemble of systems, the use of the density matrix is mandatory.

The following table lists the solvers provided by `QuantumToolbox` for dynamic quantum systems and the corresponding type of solution returned by the solver:

| **Equation** | **Function Call** | **Problem** | **Returned Solution** |
|:-------------|:------------------|:------------|:----------------------|
| Unitary evolution, Schrödinger equation | [`sesolve`](@ref) | [`sesolveProblem`](@ref) | [`TimeEvolutionSol`](@ref) |
| Lindblad master eqn. or Von Neuman eqn.| [`mesolve`](@ref) | [`mesolveProblem`](@ref) | [`TimeEvolutionSol`](@ref) |
| Monte Carlo evolution | [`mcsolve`](@ref) | [`mcsolveProblem`](@ref) [`mcsolveEnsembleProblem`](@ref) | [`TimeEvolutionMCSol`](@ref) |
| Stochastic Schrödinger equation | [`ssesolve`](@ref) | [`ssesolveProblem`](@ref) [`ssesolveEnsembleProblem`](@ref) | [`TimeEvolutionStochasticSol`](@ref) |
| Stochastic master equation | [`smesolve`](@ref) | [`smesolveProblem`](@ref) [`smesolveEnsembleProblem`](@ref) | [`TimeEvolutionStochasticSol`](@ref) |

!!! note "Solving dynamics with pre-defined problems"
    `QuantumToolbox` provides two different methods to solve the dynamics. One can use the function calls listed above by either taking all the operators (like Hamiltonian and collapse operators, etc.) as inputs directly, or generating the `prob`lems by yourself and take it as an input of the function call, e.g., `sesolve(prob)`.
