```@meta
CurrentModule = QuantumToolbox
```

# QuantumToolbox.jl Documentation

[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) is a cutting-edge Julia package designed for quantum physics simulations, closely emulating the popular Python [QuTiP](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient.

*With this package, moving from Python to Julia for quantum physics simulations has never been easier*, due to the similar syntax and functionalities.

## Features

QuantumToolbox.jl is equipped with a robust set of features:

- **Quantum State and Operator Manipulation:** Easily handle quantum states and operators with a rich set of tools, with the same functionalities as QuTiP.
- **Dynamical Evolution:** Advanced solvers for time evolution of quantum systems, thanks to the powerful [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package.
- **GPU Computing:** Leverage GPU resources for high-performance computing. For example, you run the master equation direclty on the GPU with the same syntax as the CPU case.
- **Distributed Computing:** Distribute the computation over multiple nodes (e.g., a cluster). For example, you can run undreds of quantum trajectories in parallel on a cluster, with, again, the same syntax as the simple case.
- **Easy Extension:** Easily extend the package, taking advantage of the Julia language features, like multiple dispatch and metaprogramming.

## [Installation](@id doc:Installation)

!!! note "Requirements"
    `QuantumToolbox.jl` requires `Julia 1.10+`.

To install `QuantumToolbox.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("QuantumToolbox")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.10) pkg> add QuantumToolbox
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).

To load the package and check the version information, use either [`QuantumToolbox.versioninfo()`](@ref) or [`QuantumToolbox.about()`](@ref), namely
```julia
using QuantumToolbox
QuantumToolbox.versioninfo()
QuantumToolbox.about()
```

## Brief Example

We now provide a brief example to demonstrate the similarity between [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) and [QuTiP](https://github.com/qutip/qutip).

Let's consider a quantum harmonic oscillator with a Hamiltonian given by:

```math
\hat{H} = \omega \hat{a}^\dagger \hat{a}
```

where ``\hat{a}`` and ``\hat{a}^\dagger`` are the annihilation and creation operators, respectively. We can define the Hamiltonian as follows:

```julia
using QuantumToolbox

N = 20 # cutoff of the Hilbert space dimension
ω = 1.0 # frequency of the harmonic oscillator

a = destroy(N) # annihilation operator

H = ω * a' * a
```

We now introduce some losses in a thermal environment, described by the Lindblad master equation:

```math
\frac{d \hat{\rho}}{dt} = -i [\hat{H}, \hat{\rho}] + \gamma \mathcal{D}[\hat{a}] \hat{\rho}
```

where ``\hat{\rho}`` is the density matrix, ``\gamma`` is the damping rate, and ``\mathcal{D}[\hat{a}]`` is the Lindblad dissipator, defined as:

```math
\mathcal{D}[\hat{a}]\hat{\rho} = \hat{a}\hat{\rho}\hat{a}^\dagger - \frac{1}{2}\hat{a}^\dagger\hat{a}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{a}^\dagger\hat{a}
```

We now compute the time evolution of the system using the [`mesolve`](@ref) function, starting from the initial state ``\ket{\psi (0)} = \ket{3}``:

```julia
γ = 0.1 # damping rate

ψ0 = fock(N, 3) # initial state

tlist = range(0, 10, 100) # time list

c_ops = [sqrt(γ) * a]
e_ops = [a' * a]

sol = mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops)
```

We can extract the expectation value of the number operator ``\hat{a}^\dagger \hat{a}`` with the command `sol.expect`, and the states with the command `sol.states`.

### Support for GPU calculation

We can easily pass the computation to the GPU, by simply passing all the `Qobj`s to the GPU:

```julia
using QuantumToolbox
using CUDA
CUDA.allowscalar(false) # Avoid unexpected scalar indexing

a_gpu = cu(destroy(N)) # The only difference in the code is the cu() function

H_gpu = ω * a_gpu' * a_gpu

ψ0_gpu = cu(fock(N, 3))

c_ops = [sqrt(γ) * a_gpu]
e_ops = [a_gpu' * a_gpu]

sol = mesolve(H_gpu, ψ0_gpu, tlist, c_ops, e_ops = e_ops)
```
