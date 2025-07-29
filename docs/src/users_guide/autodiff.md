# [Differentiate almost anything in QuantumToolbox](@id doc:autodiff)

Automatic differentiation (AD) has emerged as a key technique in computational science, enabling exact and efficient computation of derivatives for functions defined by code. Unlike symbolic differentiation, which may produce complex and inefficient expressions, or finite-difference methods, which suffer from numerical instability and poor scalability, AD leverages the chain rule at the level of elementary operations to provide machine-precision gradients with minimal overhead.

In `QuantumToolbox.jl`, we have introduced preliminary support for automatic differentiation. Many of the core functions are compatible with AD engines such as [`Zygote.jl`](https://github.com/FluxML/Zygote.jl), [`Enzyme.jl`](https://github.com/EnzymeAD/Enzyme.jl) or [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl), allowing users to compute gradients of observables or cost functionals involving the time evolution of open quantum systems. Although `QuantumToolbox.jl` was not originally designed with AD in mind, its architecture—rooted in Julia’s multiple dispatch and generic programming model—facilitated the integration of AD capabilities. Many core functions were already compatible with AD engines out of the box.

!!! warning "Experimental Functionality"
    At present, this functionality is considered experimental and not all parts of the library are AD-compatible. Here we provide a brief overview of the current state of AD support in `QuantumToolbox.jl` and how to use it.


::: warning

At present, this functionality is considered experimental and not all parts of the library are AD-compatible. Here we provide a brief overview of the current state of AD support in `QuantumToolbox.jl` and how to use it.

:::

## [Forward VS Reverse Mode AD](@id doc:autodiff:forward_vs_reverse)

Automatic differentiation can be broadly categorized into two modes: forward mode and reverse mode. The choice between these modes depends on the nature of the function being differentiated and the number of inputs and outputs:

- **Forward Mode AD**: This mode is particularly efficient for functions with many outputs and few inputs. It works by propagating derivatives from the inputs through the computational graph to the outputs. Forward mode is often preferred when the number of input variables is small, as it computes the derivative of each output with respect to each input in a single pass.

- **Reverse Mode AD**: In contrast, reverse mode is more efficient for functions with many inputs and few outputs. It operates by first computing the function's output and then propagating derivatives backward through the computational graph. This mode is commonly used in machine learning and optimization applications, where the loss function (output) depends on a large number of parameters (inputs).

Understanding the differences between these two modes can help users choose the most appropriate approach for their specific use case in `QuantumToolbox.jl`.

## [Differentiate the master equation](@id doc:autodiff:master_equation)

One of the primary use cases for automatic differentiation in `QuantumToolbox.jl` is the differentiation of the master equation. The master equation describes the time evolution of a quantum system's density matrix under the influence of non-unitary dynamics, such as dissipation and decoherence. Let's consider a set of parameters $\mathbf{p} = (p_1, p_2, \ldots, p_n)$ that influence the system's dynamics. The Hamiltonian and the dissipators will depend on these parameters

```math
\hat{H} = \hat{H}(\mathbf{p}), \qquad \hat{L}_j = \hat{L}_j(\mathbf{p}),
```

Hence, the density matrix will evolve according to the master equation

```math
\frac{d \hat{\rho}(\mathbf{p}, t)}{dt} = -i[\hat{H}(\mathbf{p}), \hat{\rho}(\mathbf{p}, t)] + \sum_j \hat{L}_j(\mathbf{p}) \hat{\rho}(\mathbf{p}, t) \hat{L}_j(\mathbf{p})^\dagger - \frac{1}{2} \left\{ \hat{L}_j(\mathbf{p})^\dagger \hat{L}_j(\mathbf{p}), \hat{\rho}(\mathbf{p}, t) \right\} \, ,
```

which depends on the parameters $\mathbf{p}$ and time $t$. 

We now want to compute the expectation value of an observable $\hat{O}$ at time $t$:

```math
\langle \hat{O}(\mathbf{p}, t) \rangle = \text{Tr}[\hat{O} \hat{\rho}(\mathbf{p}, t)] \, ,
```

which will also depend on the parameters $\mathbf{p}$ and time $t$.

Our goal is to compute the derivative of the expectation value with respect to the parameters:

```math
\frac{\partial \langle \hat{O}(\mathbf{p}, t) \rangle}{\partial p_j} = \frac{\partial}{\partial p_j} \text{Tr}[\hat{O} \hat{\rho}(\mathbf{p}, t)] \, .
```
