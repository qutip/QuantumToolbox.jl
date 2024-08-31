# [Quantum System Dynamics](@id doc:Quantum-System-Dynamics)

The time evolution of quantum systems lies at the heart of understanding and predicting the behavior of physical phenomena at the quantum level. Whether in theoretical research or practical applications such as quantum computing, quantum chemistry, and quantum optics, simulating the dynamics of quantum states is a fundamental task. 

Quantum systems can broadly be classified into two types: closed systems and open systems. 

- **Closed systems** are isolated from their environment and do not exchange energy or information with it. The dynamics of closed systems are unitary and reversible. The time evolution of a closed quantum system is governed by the Schrödinger equation.

- **Open systems** interact with their environment and exchange energy and information with it. The dynamics of open systems are non-unitary and irreversible. While the average dynamics of an open quantum system is described by the Lindblad master equation, a more complete description is provided by quantum trajectory approach whereby system follows a deterministic evolution conditioned on the measurement outcomes of the environment. 

The following table lists of the solvers QuantumToolbox.jl provides for dynamic quantum systems and indicates the type of solver that is returned by each function:

| Solver | Problem | Description | Return Type |
| --- | --- | --- | --- |
| [`sesolve`](@ref) | [`sesolveProblem`](@ref) | Schrödinger equation solver | [`TimeEvolutionSol`](@ref) |
| [`mesolve`](@ref) | [`mesolveProblem`](@ref) | Master equation solver | [`TimeEvolutionSol`](@ref) |
| [`lr_mesolve`](@ref) | [`lr_mesolveProblem`](@ref) | Low-rank master equation solver | [`LRTimeEvolutionSol`](@ref) |
| [`mcsolve`](@ref) | [`mcsolveProblem`](@ref) | Monte Carlo wave function solver | [`TimeEvolutionMCSol`](@ref) |
| [`ssesolve`](@ref) | [`ssesolveProblem`](@ref) | Stochastic Schrödinger equation solver | [`TimeEvolutionSSESol`](@ref) |




The following sections provide an overview of the different solvers and how to use them to simulate the time evolution of quantum systems. 

## Main differences from QuTiP

QuTip is a widely used Python library and offers grat flexibility for simulating the dynamics of quantum systems. QuantumToolbox.jl is inspired by QuTiP and aims to provide similar functionalities in Julia with close-to identical syntax. However, there are some key differences between the two libraries:

- **Performance**: QuantumToolbox.jl is built on top of the Julia programming language, which is known for its high performance and efficiency. QuantumToolbox.jl leverages Julia's just-in-time (JIT) compilation and parallel computing capabilities to vastly outperform QuTip in terms of speed and scalability. We refer to the [Performance](@ref) section for a detailed comparison of the performance of QuantumToolbox.jl and QuTip.

- **DifferentialEquations.jl**: QuantumToolbox.jl is based on DifferentialEquations.jl: a state of the art library for solving ordinary, partial, and stochastic differential equations. This library is the most efficient and flexible library for solving these problems constituting the backbone of most solvers in Julia and Python alike. The seamless integration with DifferentialEquations.jl allows QuantumToolbox.jl to take advantage of the latest advancements in numerical methods and scientific computing.

- **GPU Support**: QuantumToolbox.jl provides built-in support for GPU acceleration, allowing users to take advantage of the computational power of GPUs for simulating quantum systems.

- **Distributed Computing**: QuantumToolbox.jl supports distributed computing, enabling users with to distribute large scale problems across multiple processors or nodes for parallel execution. The structure of the solvers in QuantumToolbox.jl is designed around this feature, as discussed in the [Parallel Computing](@ref) section.

- **Callbacks**: QuantumToolbox.jl allows users to define custom callbacks that are triggered during the simulation, providing flexibility in monitoring and controlling the simulation. We encourage users to use a prefixed callback structure via the `f_ops` argument in the solvers to compute complex functions of the evolving state while avoiding the need to store the full state vector at all desired times. We refer to the [Callbacks](@ref) section for more details.

- **Additional Methods**: QuantumToolbox.jl provides additional methods and functionalities for the investigation of closed and open quantum systems. 
A few examples are: the low-rank solver [`lr_mesolve`](@ref) for the dynamics of low-entropy open quantum systems, the Arnoldi-Lindblad method [`eigsolve_al`](@ref) for faster-than-the-clock estimation of the spectrum of time-independent and Floquet open quantum systems, the breath-first seach algorithm [`bdf`](@ref) for the block-diagonalization of aribitrary Hamiltonians/Liouvillians with unknown symmetries, and the U(1) symmetric Monte-Carlo solver [`mcsolve`](@ref).