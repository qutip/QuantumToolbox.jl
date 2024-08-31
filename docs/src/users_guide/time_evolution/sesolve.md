# Unitary dynamics

The evolution of a closed quantum system is governed by the Schrödinger equation
```math
i \frac{d}{dt} \left| \psi(t) \right\rangle = \hat{H} \left| \psi(t) \right\rangle,
```
where $\hat{H}$ is the Hamiltonian of the system, and $\left| \psi(t) \right\rangle$ is the state vector. 
Numerically, $\left| \psi(t) \right\rangle$ is represented as a column vector of size $N\times 1$, while $\hat{H}$ as a matrix of size $N\times N$. $N$ is the dimension of the Hilbert space scaling exponentially with the number of degrees of freedom in the system (e.g., the number of qubits in a quantum circuit). 
 The formal solution to the Schrödinger equation is given by
```math 
|\psi(t)\rangle = e^{-i\hat{H}t}|\psi(0)\rangle.
```
The Hermiticity of $\hat{H}$, makes the evolution operation unitary and the dynamics reversible. The norm of the state vector, as well as the expectation value of $\hat{H}$ are preserved under unitary evolution.

Analytic solutions to the Schrödinger equation are rare and limited to very small systems making numerical methods the only practical approach for larger systems. A direct approach is to directly compute the matrix exponential of $\hat{H}$. This scales as $O(N^3)$, similar to matrix diagonalization, making it impractical for large systems. In contrast, matrix-vector multiplication scales as $O(N^2)$.
Integration schemes for the Schrödinger equation involve solving the equation directly using matrix-vector products with adaptive time steps and error control. The Euler method is a simple example of a numerical integration scheme that approximates the time evolution over an infinitesimal time step $\Delta t$ as
```math
|\psi(t+\Delta t)\rangle = (1-i\hat{H}\Delta t)|\psi(t)\rangle
```
More sophisticated schemes like Runge-Kutta or Crank-Nicolson methods are often used in practice to achieve better accuracy and stability.

QuantumToolbox.jl provides the `sesolve` function to simulate the unitary dynamics of a quantum system. The function takes a `sesolveProblem` object as input and returns a `TimeEvolutionSol` object that contains the solution to the Schrödinger equation. The following example demonstrates how to use the `sesolve` function to simulate the time evolution of a simple quantum system.