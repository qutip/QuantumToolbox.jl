# [Hierarchical Equations of Motion](@id doc:Hierarchical-Equations-of-Motion)

The hierarchical equations of motion (HEOM) approach was originally developed by [Tanimura1989](@citet) in the context of physical chemistry to "exactly" solve a quantum system (labeled as ``\textrm{s}``) in contact with a bosonic environment, encapsulated in the following total Hamiltonian:

```math
\hat{H}_{\textrm{total}} = \hat{H}_{\textrm{s}} + \sum_k \omega_k \hat{b}^\dagger_k \hat{b}_k + \hat{V}_{\textrm{s}} \sum_k g_k \left(\hat{b}_k + \hat{b}^\dagger_k\right),
```

where ``\hat{b}_k`` (``\hat{b}^\dagger_k``) is the bosonic annihilation (creation) operator associated to the ``k``th mode (with frequency ``\omega_k``), ``\hat{V}_{\textrm{s}}`` refer to the coupling operator acting on the system's degree of freedom, and ``g_k`` are the coupling strengths.

As in other solutions to this problem, the properties of the bath are encapsulated by its temperature and its spectral density,

```math
J(\omega) = 2 \pi \sum_k g^2_k \delta(\omega - \omega_k).
```

In the HEOM approach, for bosonic baths, one typically chooses a Drude-Lorentz spectral density:

```math
J_{\textrm{DL}}(\omega) = \frac{4 \Delta W \omega}{\omega^2 + W^2},
```

or an under-damped Brownian motion spectral density,

```math
J_{\textrm{U}}(\omega)=\frac{2 \Delta^2 W \omega}{(\omega^2 - \omega_0^2)^2 + \omega^2 W^2}.
```

Here, ``\Delta`` represents the coupling strength between the system and the bosonic bath with band-width ``W`` and resonance frequency ``\omega_0``.

We introduce an efficient `Julia` framework for HEOM approach called [`HierarchicalEOM.jl`](https://github.com/qutip/HierarchicalEOM.jl). This package is built upon `QuantumToolbox.jl` and provides a user-friendly and efficient tool to simulate complex open quantum systems based on HEOM approach. For a detailed explanation of this package, we recommend to read its [documentation](https://qutip.org/HierarchicalEOM.jl/) and also the article [Huang2023](@citet).

Given the spectral density, the HEOM approach requires a decomposition of the bath correlation functions in terms of exponentials. In the [documentation of `HierarchicalEOM.jl`](https://qutip.org/HierarchicalEOM.jl/), we not only describe how this is done for both bosonic and fermionic environments with code examples, but also describe how to solve the time evolution (dynamics), steady-states, and spectra based on HEOM approach.
