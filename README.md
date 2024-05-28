<!-- Logo -->
<div align="center">
  <img src="./docs/src/assets/logo.png" alt="QuantumToolbox.jl logo" width="150">
</div>

# QuantumToolbox.jl

<!-- Current admin team (in chronological order) -->
[A. Mercurio](https://github.com/albertomercurio),
[L. Gravina](https://github.com/lgravina1997),
and [Y.-T. Huang](https://github.com/ytdHuang).

<!-- Badges Table -->
| **Release**       | [![Release][release-img]](release-url) [![License][license-img]](license-url) [![DOI][doi-img]](doi-url) [![Downloads][download-img]](download-url) |
|:-----------------:|:-------------|
| **Runtests**      | [![Runtests][runtests-img]](runtests-url) [![Coverage][codecov-img]](codecov-url) [![Aqua QA][aqua-img]](aqua-url) [![JET][jet-img]](jet-url) |
| **Documentation** | [![Doc-Stable][docs-stable-img]](docs-stable-url) [![Doc-Dev][docs-develop-img]](docs-develop-url) |

[release-img]: https://img.shields.io/github/release/qutip/QuantumToolbox.jl.svg
[release-url]: https://github.com/qutip/QuantumToolbox.jl/releases

[license-img]: https://img.shields.io/badge/license-New%20BSD-blue.svg
[license-url]: https://opensource.org/licenses/BSD-3-Clause

[doi-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.10822816.svg
[doi-url]: https://doi.org/10.5281/zenodo.10822816

[download-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FQuantumToolbox&query=total_requests&label=Downloads
[download-url]: https://juliapkgstats.com/pkg/QuantumToolbox

[runtests-img]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/CI.yml/badge.svg?branch=main
[runtests-url]: https://github.com/qutip/QuantumToolbox.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/qutip/QuantumToolbox.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/qutip/QuantumToolbox.jl

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[jet-img]: https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-9cf
[jet-url]: https://github.com/aviatesk/JET.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://qutip.github.io/QuantumToolbox.jl/stable
[docs-develop-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-develop-url]: https://qutip.github.io/QuantumToolbox.jl/dev

## Introduction

[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) is a cutting-edge Julia package designed for quantum physics simulations, closely emulating the popular Python [QuTiP](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient.

## Features

QuantumToolbox.jl is equipped with a robust set of features:

- **Quantum State and Operator Manipulation:** Easily handle quantum states and operators with a rich set of tools.
- **Dynamical Evolution:** Advanced solvers for time evolution of quantum systems.
- **Measurement and Statistics:** Comprehensive quantum measurement simulation and analysis.
- **GPU and Distributed Computing:** Leverage GPU and distributed resources for high-performance computing.

## Installation
`QuantumToolbox.jl` requires `Julia 1.7+`. To install it, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("QuantumToolbox")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.7) pkg> add QuantumToolbox
```
To load the package and check the version information, use either `versioninfo()` or `about()`, namely
```julia
using QuantumToolbox
QuantumToolbox.versioninfo()
QuantumToolbox.about()
```
