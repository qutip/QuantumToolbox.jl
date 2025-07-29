# [Automatic Differentiation](@id doc:autodiff)  

Automatic differentiation (AD) has emerged as a key technique in computational science, enabling exact and efficient computation of derivatives for functions defined by code. Unlike symbolic differentiation, which may produce complex and inefficient expressions, or finite-difference methods, which suffer from numerical instability and poor scalability, AD leverages the chain rule at the level of elementary operations to provide machine-precision gradients with minimal overhead.

In `QuantumToolbox.jl`, we have introduced preliminary support for automatic differentiation. Many of the core functions are compatible with AD engines such as [`Zygote.jl`](https://github.com/FluxML/Zygote.jl), [`Enzyme.jl`](https://github.com/EnzymeAD/Enzyme.jl) or [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl), allowing users to compute gradients of observables or cost functionals involving the time evolution of open quantum systems. Although `QuantumToolbox.jl` was not originally designed with AD in mind, its architecture—rooted in Julia’s multiple dispatch and generic programming model—facilitated the integration of AD capabilities. Many core functions were already compatible with AD engines out of the box.

!!! warning "Experimental Functionality"
    At present, this functionality is considered experimental and not all parts of the library are AD-compatible. Here we provide a brief overview of the current state of AD support in `QuantumToolbox.jl` and how to use it.


## [Forward versus Reverse Mode AD](@id doc:autodiff:forward-versus-reverse) 

Automatic differentiation can be broadly categorized into two modes: forward mode and reverse mode. The choice between these modes depends on the nature of the function being differentiated and the number of inputs and outputs:

- **Forward Mode AD**: This mode is particularly efficient for functions with many outputs and few inputs. It works by propagating derivatives from the inputs through the computational graph to the outputs. Forward mode is often preferred when the number of input variables is small, as it computes the derivative of each output with respect to each input in a single pass.

- **Reverse Mode AD**: In contrast, reverse mode is more efficient for functions with many inputs and few outputs. It operates by first computing the function's output and then propagating derivatives backward through the computational graph. This mode is commonly used in machine learning and optimization applications, where the loss function (output) depends on a large number of parameters (inputs).

Understanding the differences between these two modes can help users choose the most appropriate approach for their specific use case in `QuantumToolbox.jl`.

## [Differentiate the master equation](@id doc:autodiff:master-equation)

One of the primary use cases for automatic differentiation in `QuantumToolbox.jl` is the differentiation of the master equation. The master equation describes the time evolution of a quantum system's density matrix under the influence of non-unitary dynamics, such as dissipation and decoherence. Let's consider a set of parameters $\mathbf{p} = (p_1, p_2, \ldots, p_n)$ that influence the system's dynamics. The Hamiltonian and the dissipators will depend on these parameters

```math
\hat{H} = \hat{H}(\mathbf{p}), \qquad \hat{L}_j = \hat{L}_j(\mathbf{p}),
```

Hence, the density matrix will evolve according to the master equation

```@raw html
<span id="eq:master-equation"></span>
```
```math
\begin{align}
\frac{d \hat{\rho}(\mathbf{p}, t)}{dt} =& -i[\hat{H}(\mathbf{p}), \hat{\rho}(\mathbf{p}, t)] \\
&+ \sum_j \hat{L}_j(\mathbf{p}) \hat{\rho}(\mathbf{p}, t) \hat{L}_j(\mathbf{p})^\dagger - \frac{1}{2} \left\{ \hat{L}_j(\mathbf{p})^\dagger \hat{L}_j(\mathbf{p}), \hat{\rho}(\mathbf{p}, t) \right\} \, , 
\end{align} \tag{1}
```

which depends on the parameters $\mathbf{p}$ and time $t$. 

We now want to compute the expectation value of an observable $\hat{O}$ at time $t$:

```math
\langle \hat{O}(\mathbf{p}, t) \rangle = \text{Tr}[\hat{O} \hat{\rho}(\mathbf{p}, t)] \, ,
```

which will also depend on the parameters $\mathbf{p}$ and time $t$.

Our goal is to compute the derivative of the expectation value with respect to the parameters:

```math
\frac{\partial \langle \hat{O}(\mathbf{p}, t) \rangle}{\partial p_j} = \frac{\partial}{\partial p_j} \text{Tr}[\hat{O} \hat{\rho}(\mathbf{p}, t)] \, ,
```

and to achieve this, we can use an AD engine like [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) (forward mode) or [`Zygote.jl`](https://github.com/FluxML/Zygote.jl) (reverse mode).

Let's apply this to a simple example of a driven-dissipative quantum harmonic oscillator. The Hamiltonian in the drive frame is given by

```math
\hat{H} = \Delta \hat{a}^\dagger \hat{a} + F \left( \hat{a} + \hat{a}^\dagger \right) \, ,
```

where $\Delta = \omega_0 - \omega_d$ is the cavity-drive detuning, $F$ is the drive strength, and $\hat{a}$ and $\hat{a}^\dagger$ are the annihilation and creation operators, respectively. The system is subject to a single dissipative channel with a Lindblad operator $\hat{L} = \sqrt{\gamma} \hat{a}$, where $\gamma$ is the dissipation rate. If we start from the ground state $\hat{\rho}(0) = \vert 0 \rangle \langle 0 \vert$, the systems evolves according to the master equation in [Eq. (1)](#eq:master-equation).

We now want to study the number of photons at the steady state, and how it varies with $\mathbf{p} = (\Delta, F, \gamma)$, namely $\nabla_\mathbf{p} \langle \hat{a}^\dagger \hat{a} \rangle (\mathbf{p}, t \to \infty)$. We can extract an analytical expression, in order to verify the correctness of the AD implementation:

```math
\langle \hat{a}^\dagger \hat{a} \rangle_\mathrm{ss} = \frac{F^2}{\Delta^2 + \frac{\gamma^2}{4}} \, ,
```

with the gradient given by

```math
\nabla_\mathbf{p} \langle \hat{a}^\dagger \hat{a} \rangle_\mathrm{ss} =
\begin{pmatrix}
\frac{-2 F^2 \Delta}{(\Delta^2 + \frac{\gamma^2}{4})^2} \\
\frac{2 F}{\Delta^2 + \frac{\gamma^2}{4}} \\
\frac{-F^2 \gamma}{2 (\Delta^2 + \frac{\gamma^2}{4})^2}
\end{pmatrix} \, .
```

Although `QuantumToolbox.jl` has the [`steadystate`](@ref) function to directly compute the steady state without explicitly solving the master equation, here we use the [`mesolve`](@ref) function to integrate up to a long time $t_\mathrm{max}$, and then compute the expectation value of the number operator. We will demonstrate how to compute the gradient using both [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and [`Zygote.jl`](https://github.com/FluxML/Zygote.jl).

### [Forward Mode AD with ForwardDiff.jl](@id doc:autodiff:forward)

```@setup autodiff
using QuantumToolbox
```

We start by importing [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and defining the parameters and operators:

```@example autodiff
using ForwardDiff

const N = 20
const a = destroy(N)
const ψ0 = fock(N, 0)
const t_max = 40
const tlist = range(0, t_max, 100)
```

Then, we define a function that take the parameters `p` as an input and returns the expectation value of the number operator at `t_max`. We also define the analytical solution of the steady state photon number and its gradient for comparison:

```@example autodiff
function my_f_mesolve_direct(p)
    H = p[1] * a' * a + p[2] * (a + a')
    c_ops = [sqrt(p[3]) * a]
    sol = mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
    return real(expect(a' * a, sol.states[end]))
end

# Analytical solution
function my_f_analytical(p)
    Δ, F, γ = p
    return F^2 / (Δ^2 + γ^2 / 4)
end
function my_grad_analytical(p)
    Δ, F, γ = p
    return [
        -2 * F^2 * Δ / (Δ^2 + γ^2 / 4)^2,
         2 * F / (Δ^2 + γ^2 / 4),
        -F^2 * γ / (2 * (Δ^2 + γ^2 / 4)^2)
    ]
end
```

The gradient can be computed using `ForwardDiff.gradient`:

```@example autodiff
Δ = 1.5
F = 1.5
γ = 1.5
params = [Δ, F, γ]

grad_exact = my_grad_analytical(params)
grad_fd = ForwardDiff.gradient(my_f_mesolve_direct, params)
```

and test if the results match:

```@example autodiff
isapprox(grad_exact, grad_fd; atol = 1e-5)
```

### [Reverse Mode AD with Zygote.jl](@id doc:autodiff:reverse)

Reverse-mode differentiation is significantly more challenging than forward-mode when dealing ODEs, as the complexity arises from the need to propagate gradients backward through the entire time evolution of the quantum state.

`QuantumToolbox.jl` leverages the advanced capabilities of [`SciMLSensitivity.jl`](https://github.com/SciML/SciMLSensitivity.jl) to handle this complexity. [`SciMLSensitivity.jl`](https://github.com/SciML/SciMLSensitivity.jl) implements sophisticated methods for computing gradients of ODE solutions, such as the adjoint method, which computes gradients by solving an additional "adjoint" ODE backward in time. For more details on the adjoint method and other sensitivity analysis techniques, please refer to the [`SciMLSensitivity.jl` documentation](https://docs.sciml.ai/SciMLSensitivity/stable/).

In order to reverse-differentiate the master equation, we need to define the operators as [`QuantumObjectEvolution`](@ref) objects, which use [`SciMLOperators.jl`](https://github.com/SciML/SciMLOperators.jl) to represent parameter-dependent operators.

```@example autodiff
using Zygote
using SciMLSensitivity

# For SciMLSensitivity.jl
coef_Δ(p, t) = p[1]
coef_F(p, t) = p[2]
coef_γ(p, t) = sqrt(p[3])
H = QobjEvo(a' * a, coef_Δ) + QobjEvo(a + a', coef_F)
c_ops = [QobjEvo(a, coef_γ)]
const L = liouvillian(H, c_ops)

function my_f_mesolve(p)
    sol = mesolve(
        L,
        ψ0,
        tlist,
        progress_bar = Val(false),
        params = p,
        sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()),
    )

    return real(expect(a' * a, sol.states[end]))
end
```

And the gradient can be computed using `Zygote.gradient`:

```@example autodiff
grad_zygote = Zygote.gradient(my_f_mesolve, params)[1]
```

Finally, we can compare the results from [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and [`Zygote.jl`](https://github.com/FluxML/Zygote.jl):

```@example autodiff
isapprox(grad_fd, grad_zygote; atol = 1e-5)
```

## [Conclusion](@id doc:autodiff:conclusion)

In this section, we have explored the integration of automatic differentiation into `QuantumToolbox.jl`, enabling users to compute gradients of observables and cost functionals involving the time evolution of open quantum systems. We demonstrated how to differentiate the master equation using both forward mode with [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and reverse mode with [`Zygote.jl`](https://github.com/FluxML/Zygote.jl), showcasing the flexibility and power of automatic differentiation in quantum computing applications. AD can be applied to other functions in `QuantumToolbox.jl`, although the support is still experimental and not all functions are guaranteed to be compatible. We encourage users to experiment with AD in their quantum simulations and contribute to the ongoing development of this feature.