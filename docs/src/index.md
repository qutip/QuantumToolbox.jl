```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "QuantumToolbox.jl"
  tagline: A pure Julia framework designed for High-performance quantum physics simulations
  image:
    src: /logo.png
    alt: QuantumToolbox
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: Users Guide
      link: /users_guide/QuantumObject/QuantumObject
    - theme: alt
      text: View on Github
      link: https://github.com/qutip/QuantumToolbox.jl
    - theme: alt
      text: API
      link: /resources/api
    - theme: alt
      text: Visit QuTiP.org
      link: https://qutip.org/


features:
  - icon: <img width="64" height="64" src="https://docs.sciml.ai/DiffEqDocs/stable/assets/logo.png" alt="markdown"/>
    title: Dynamical Evolution
    details: Advanced solvers for time evolution of quantum systems, thanks to the powerful DifferentialEquations.jl package.
    link: /users_guide/time_evolution/intro
  - icon: <img width="64" height="64" src="https://cuda.juliagpu.org/stable/assets/logo.png" />
    title: GPU Computing
    details: Leverage GPU resources for high-performance computing. Simulate quantum dynamics directly on the GPU with the same syntax as the CPU case.
    link: /getting_started
  - icon: <img width="64" height="64" src="https://img.icons8.com/?size=100&id=1W4Bkj363ov0&format=png&color=000000" />
    title: Distributed Computing
    details: Distribute the computation over multiple nodes (e.g., a cluster). Simulate hundreds of quantum trajectories in parallel on a cluster, with, again, the same syntax as the simple case.
    link: /users_guide/time_evolution/mcsolve
---
```

# [Introduction](@id doc:Introduction)

[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) is a cutting-edge Julia package designed for quantum physics simulations, closely emulating the popular Python [QuTiP](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient. Taking advantage of the Julia language features (like multiple dispatch and metaprogramming), QuantumToolbox.jl is designed to be easily extendable, allowing users to build upon the existing functionalities.

*__With this package, moving from Python to Julia for quantum physics simulations has never been easier__*, due to the similar syntax and functionalities.

# [Installation](@id doc:Installation)

!!! note "Requirements"
    `QuantumToolbox.jl` requires `Julia 1.10+`.

To install `QuantumToolbox.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("QuantumToolbox")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-repl
(1.10) pkg> add QuantumToolbox
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).

To load the package and check the version information, use either [`QuantumToolbox.versioninfo()`](@ref) or [`QuantumToolbox.about()`](@ref), namely
```julia
using QuantumToolbox
QuantumToolbox.versioninfo()
QuantumToolbox.about()
```
