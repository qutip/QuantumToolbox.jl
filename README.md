<!-- Logo -->
<div align="center">
  <img src="./docs/src/assets/logo.png" alt="QuantumToolbox.jl logo" width="150">
</div>

# QuantumToolbox.jl

<!-- Current admin team (in chronological order) -->
[A. Mercurio](https://github.com/albertomercurio)
and [Y.-T. Huang](https://github.com/ytdHuang).

<!-- Table of Badges -->
| **Release**       | [![Release][release-img]][release-url] [![License][license-img]][license-url] [![Cite][cite-img]][cite-url] [![Downloads][download-img]][download-url] |
|:-----------------:|:-------------|
| **Runtests**      | [![Runtests][runtests-img]][runtests-url] [![Coverage][codecov-img]][codecov-url] |
| **Code Quality**  | [![Code Quality][code-quality-img]][code-quality-url] [![Aqua QA][aqua-img]][aqua-url] [![JET][jet-img]][jet-url] |
| **Documentation** | [![Doc-Stable][docs-stable-img]][docs-stable-url] [![Doc-Dev][docs-develop-img]][docs-develop-url] |
| **Benchmark** | [![Benchmarks][benchmark-img]][benchmark-url] |
| **Community** | [![Zulip][zulip-img]][zulip-url] [![QuTiP-discussion][QuTiP-discussion-img]][QuTiP-discussion-url] |
| **Support** | [![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](https://unitary.fund) |

[release-img]: https://img.shields.io/github/release/qutip/QuantumToolbox.jl.svg
[release-url]: https://github.com/qutip/QuantumToolbox.jl/releases

[license-img]: https://img.shields.io/badge/license-New%20BSD-blue.svg
[license-url]: https://opensource.org/licenses/BSD-3-Clause

[cite-img]: https://img.shields.io/badge/cite-Quantum_9%2C_1866_(2025)-blue
[cite-url]: https://doi.org/10.22331/q-2025-09-29-1866

[download-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FQuantumToolbox&query=total_requests&label=Downloads
[download-url]: https://juliapkgstats.com/pkg/QuantumToolbox

[runtests-img]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/CI.yml/badge.svg?branch=main
[runtests-url]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/qutip/QuantumToolbox.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/qutip/QuantumToolbox.jl

[code-quality-img]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/Code-Quality.yml/badge.svg 
[code-quality-url]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/Code-Quality.yml

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[jet-img]: https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a
[jet-url]: https://github.com/aviatesk/JET.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://qutip.github.io/QuantumToolbox.jl/stable
[docs-develop-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-develop-url]: https://qutip.github.io/QuantumToolbox.jl/dev

[benchmark-img]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/Benchmarks.yml/badge.svg?branch=main
[benchmark-url]: https://qutip.org/QuantumToolbox.jl/benchmarks/

[zulip-img]: https://img.shields.io/badge/Zulip%20Chat-join-6f73af.svg
[zulip-url]: https://quantumtoolbox-jl.zulipchat.com

[QuTip-discussion-img]: https://img.shields.io/badge/QuTiP%20discussion%20group-join-6f73af.svg
[QuTip-discussion-url]: https://groups.google.com/g/qutip

## Introduction

[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) is a cutting-edge [`Julia`](https://julialang.org/) package designed for quantum physics simulations, closely emulating the popular Python [`QuTiP`](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of [`Julia`](https://julialang.org/) with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient.

*With this package, moving from Python to Julia for quantum physics simulations has never been easier*, due to the similar syntax and functionalities.

## Features

`QuantumToolbox.jl` is equipped with a robust set of features:

- **Quantum State and Operator Manipulation:** Easily handle quantum states and operators with a rich set of tools, with the same functionalities as `QuTiP`.
- **Dynamical Evolution:** Advanced solvers for time evolution of quantum systems, thanks to the powerful [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) package.
- **GPU Computing:** Leverage GPU resources for high-performance computing. Simulate quantum dynamics directly on the GPU with the same syntax as the CPU case.
- **Distributed Computing:** Distribute the computation over multiple nodes (e.g., a cluster). For example, you can run hundreds of quantum trajectories in parallel on a cluster, with, again, the same syntax as the simple case. See [here](https://qutip.org/QuantumToolbox.jl/stable/users_guide/cluster) for more information.
- **Differentiable Programming:** Enable gradient-based optimization for quantum algorithms. Compute gradients of quantum dynamics with respect to their parameters using automatic differentiation. See [here](https://qutip.org/QuantumToolbox.jl/stable/users_guide/autodiff) for more information.
- **Easy Extension:** Easily extend the package, taking advantage of the `Julia` language features, like multiple dispatch and metaprogramming.

## Installation
    
> [!NOTE]
> `QuantumToolbox.jl` requires `Julia 1.10+`.

To install `QuantumToolbox.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("QuantumToolbox")
```
Alternatively, this can also be done in `Julia`'s [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-repl
(1.10) pkg> add QuantumToolbox
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).

To load the package and check the version information, use either `QuantumToolbox.versioninfo()` or `QuantumToolbox.about()`, namely
```julia
using QuantumToolbox
QuantumToolbox.versioninfo()
QuantumToolbox.about()
```

## Brief Example

We now provide a brief example to demonstrate the similarity between [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) and [`QuTiP`](https://github.com/qutip/qutip).

Let's consider a quantum harmonic oscillator with a Hamiltonian given by:

$$
\hat{H} = \omega \hat{a}^\dagger \hat{a}
$$

where $\hat{a}$ and $\hat{a}^\dagger$ are the annihilation and creation operators, respectively. We can define the Hamiltonian as follows:

```julia
using QuantumToolbox

N = 20 # cutoff of the Hilbert space dimension
ω = 1.0 # frequency of the harmonic oscillator

a = destroy(N) # annihilation operator

H = ω * a' * a
```

We now introduce some losses in a thermal environment, described by the Lindblad master equation:

$$
\frac{d \hat{\rho}}{dt} = -i [\hat{H}, \hat{\rho}] + \gamma \mathcal{D}[\hat{a}] \hat{\rho}
$$

where $\hat{\rho}$ is the density matrix, $\gamma$ is the damping rate, and $\mathcal{D}[\hat{a}]$ is the Lindblad dissipator, defined as:

$$
\mathcal{D}[\hat{a}]\hat{\rho} = \hat{a}\hat{\rho}\hat{a}^\dagger - \frac{1}{2}\hat{a}^\dagger\hat{a}\hat{\rho} - \frac{1}{2}\hat{\rho}\hat{a}^\dagger\hat{a}
$$

We now compute the time evolution of the system using the `mesolve` function, starting from the initial state $\ket{\psi (0)} = \ket{3}$:

```julia
γ = 0.1 # damping rate

ψ0 = fock(N, 3) # initial state

tlist = range(0, 10, 100) # time list

c_ops = [sqrt(γ) * a]
e_ops = [a' * a]

sol = mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops)
```

We can extract the expectation value of the number operator $\hat{a}^\dagger \hat{a}$ with the command `sol.expect`, and the states with the command `sol.states`.

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

## Performance comparison with other packages

Here we provide a brief performance comparison between `QuantumToolbox.jl` and other popular quantum physics simulation packages, such as [`QuTiP`](https://github.com/qutip/qutip) (Python), [`dynamiqs`](https://github.com/dynamiqs/dynamiqs) (Python - JAX) and [`QuantumOptics.jl`](https://github.com/qojulia/QuantumOptics.jl) (Julia). We clearly show that `QuantumToolbox.jl` is the fastest package among the four. A detailed code is available [here](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/benchmarks.jl).

![](https://raw.githubusercontent.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/refs/heads/main/figures/benchmarks.svg)

## Contributing to QuantumToolbox.jl

You are most welcome to contribute to `QuantumToolbox.jl` development by forking this repository and sending pull requests (PRs), or filing bug reports at the issues page.

Contributors and users for `QuantumToolbox.jl` are invited to [![Zulip][zulip-img]][zulip-url] for questions, development discussions, and community updates.

You can also help out with users' questions, or discuss proposed changes for the entire QuTiP organization in the [![QuTiP-discussion][QuTiP-discussion-img]][QuTiP-discussion-url].

For more information about contribution, including technical advice, please see the [Contributing to Quantum Toolbox in Julia](https://qutip.org/QuantumToolbox.jl/stable/resources/contributing).

## Cite `QuantumToolbox.jl`
If you like `QuantumToolbox.jl`, we would appreciate it if you starred the repository in order to help us increase its visibility. Furthermore, if you find the framework useful in your research, we would be grateful if you could cite our publication [ [Quantum 9, 1866 (2025)](https://doi.org/10.22331/q-2025-09-29-1866)  ] using the following bibtex entry:

```bib
@article{QuantumToolbox.jl2025,
  title = {Quantum{T}oolbox.jl: {A}n efficient {J}ulia framework for simulating open quantum systems},
  author = {Mercurio, Alberto and Huang, Yi-Te and Cai, Li-Xun and Chen, Yueh-Nan and Savona, Vincenzo and Nori, Franco},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {9},
  pages = {1866},
  month = sep,
  year = {2025},
  doi = {10.22331/q-2025-09-29-1866},
  url = {https://doi.org/10.22331/q-2025-09-29-1866}
}
```

## Acknowledgements

### Fundings

`QuantumToolbox.jl` is supported by the [Unitary Fund](https://unitary.fund), a grant program for quantum technology projects.

<div align="center">
  <a href="https://unitary.fund" target="about:blank">
    <img src="https://raw.githubusercontent.com/unitaryfund/unitary.fund/refs/heads/main/src/assets/svg/logo.svg" alt="Unitary Fund logo" width="200">
  </a>
</div>

`QuantumToolbox.jl` is partially supported by

<div align="center">
  <a href="https://www.epfl.ch" target="about:blank"><img src="docs/src/assets/org-logo/EPFL.svg" style="margin:5px;" alt="EPFL logo" width="200"></a> 
  <a href="https://www.epfl.ch/research/domains/quantum-center" target="about:blank"><img src="docs/src/assets/org-logo/EPFL-QSE.svg" style="margin:5px;" alt="EPFL QSE logo" width="200"></a> 

  <a href="https://www.ncku.edu.tw" target="about:blank"><img src="docs/src/assets/org-logo/NCKU.png" style="margin:5px;" alt="NCKU logo" width="250"></a> 
  <a href="https://qfort.ncku.edu.tw" target="about:blank"><img src="docs/src/assets/org-logo/QFort.png" style="margin:5px;" alt="QFort logo" width="250"></a> 
  <a href="https://www.nstc.gov.tw" target="about:blank"><img src="docs/src/assets/org-logo/NSTC.png" style="margin:5px;" alt="NSTC logo" width="250"></a> 

  <a href="https://dml.riken.jp" target="about:blank"><img src="docs/src/assets/org-logo/RIKEN.png" style="margin:5px;" alt="RIKEN logo" width="200"></a> 
  <a href="https://rqc.riken.jp" target="about:blank"><img src="docs/src/assets/org-logo/RIKEN-RQC.svg" style="margin:5px;" alt="RIKEN RQC logo" width="200"></a>
  <a href="https://www.jst.go.jp/moonshot" target="about:blank"><img src="docs/src/assets/org-logo/JST-Moonshot.png" style="margin:5px;" alt="JST Moonshot logo" width="200"></a>
</div>

### Other Acknowledgements

We are also grateful to the [Zulip](https://zulip.com) team for providing a free chat service for open-source projects.

<div align="center">
  <a href="https://zulip.com" target="about:blank">
    <img src="https://zulip.com/static/images/logo/zulip-org-logo.svg" alt="Zulip logo" width="200">
  </a>
</div>
