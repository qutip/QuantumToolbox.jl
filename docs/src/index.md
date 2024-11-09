```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "QuantumToolbox.jl"
  tagline: High-performance quantum simulations made simple
  image:
    src: /logo.png
    alt: QuantumToolbox
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/qutip/QuantumToolbox.jl
    - theme: alt
      text: API
      link: /api


features:
  - icon: <img width="64" height="64" src="https://docs.sciml.ai/DiffEqDocs/stable/assets/logo.png" alt="markdown"/>
    title: Dynamical Evolution
    details: Advanced solvers for time evolution of quantum systems, thanks to the powerful DifferentialEquations.jl package.
    link: /users_guide/time_evolution/intro
  - icon: <img width="64" height="64" src="https://cuda.juliagpu.org/stable/assets/logo.png" />
    title: GPU Computing
    details: Leverage GPU resources for high-performance computing. Simulate the master equation directly on the GPU with the same syntax as the CPU case.
    link: /getting_started
  - icon: <img width="64" height="64" src="https://img.icons8.com/?size=100&id=1W4Bkj363ov0&format=png&color=000000" />
    title: Distributed Computing
    details: Distribute the computation over multiple nodes (e.g., a cluster). Simulate undreds of quantum trajectories in parallel on a cluster, with, again, the same syntax as the simple case.
    link: /users_guide/time_evolution/mcsolve
---
```

[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) is a cutting-edge Julia package designed for quantum physics simulations, closely emulating the popular Python [QuTiP](https://github.com/qutip/qutip) package. It uniquely combines the simplicity and power of Julia with advanced features like GPU acceleration and distributed computing, making simulation of quantum systems more accessible and efficient. Taking advantage of the Julia language features (like multiple dispatch and metaprogramming), QuantumToolbox.jl is designed to be easily extendable, allowing users to build upon the existing functionalities.

*__With this package, moving from Python to Julia for quantum physics simulations has never been easier__*, due to the similar syntax and functionalities.
