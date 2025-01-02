```@meta
CurrentModule = QuantumToolbox
```

# Brief Example

We now provide a brief example to demonstrate the similarity between [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) and [`QuTiP`](https://github.com/qutip/qutip).

## CPU Computation

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

!!! note "Lindblad master equation"
    See [here](@ref doc-TE:Lindblad-Master-Equation-Solver) for more details about Lindblad master equation.

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

## GPU Computation

!!! note "Extension for CUDA.jl"
    `QuantumToolbox.jl` provides an extension to support GPU computation. To trigger the extension, you need to install and import [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl) together with `QuantumToolbox.jl`. See [here](@ref doc:CUDA) for more details.

```julia
using QuantumToolbox
using CUDA
CUDA.allowscalar(false) # Avoid unexpected scalar indexing
```

We can easily pass the computation to the GPU, by simply passing all the [`QuantumObject`](@ref)s to the GPU:

```julia
a_gpu = cu(destroy(N)) # The only difference in the code is the cu() function

H_gpu = ω * a_gpu' * a_gpu

ψ0_gpu = cu(fock(N, 3))

c_ops = [sqrt(γ) * a_gpu]
e_ops = [a_gpu' * a_gpu]

sol = mesolve(H_gpu, ψ0_gpu, tlist, c_ops, e_ops = e_ops)
```