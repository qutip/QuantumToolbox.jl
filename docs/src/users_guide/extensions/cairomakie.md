# [Extension for the Makie.jl ecosystem](@id doc:Makie)

This is an extension to support visualization (plotting functions) using [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) library.

This extension will be automatically loaded if user imports both `QuantumToolbox.jl` and [`Makie.jl`](https://github.com/MakieOrg/Makie.jl). It is worth noting that the `Makie.jl` package provides only the engine for plotting, and the user has to import the specific backend. Here we demonstrate the usage of [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie), which will automatically import `Makie.jl`.

```julia
using QuantumToolbox
using CairoMakie
```

To plot with [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie) library, specify the keyword argument `library = Val(:Makie)` for the plotting functions.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:Makie)` instead of `:Makie`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.

The supported plotting functions are listed as follows:

| **Plotting Function** | **Description** |
|:----------------------|:----------------|
| [`plot_wigner`](@ref) | [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) |
| [`plot_fock_distribution`](@ref) | [Fock state](https://en.wikipedia.org/wiki/Fock_state) distribution |
| [`plot_bloch`](@ref)  | [Plotting on the Bloch Sphere](https://qutip.readthedocs.io/en/stable/guide/guide-bloch.html) |
