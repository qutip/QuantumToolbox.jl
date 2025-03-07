# [Extension for CairoMakie.jl](@id doc:CairoMakie)

This is an extension to support visualization (plotting functions) using [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie) library.

This extension will be automatically loaded if user imports both `QuantumToolbox.jl` and [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie):

```julia
using QuantumToolbox
using CairoMakie
```

To plot with [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie) library, specify the keyword argument `library = Val(:CairoMakie)` for the plotting functions.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:CairoMakie)` instead of `:CairoMakie`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.

The supported plotting functions are listed as follows:

| **Plotting Function** | **Description** |
|:----------------------|:----------------|
| [`plot_wigner`](@ref) | [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) |
| [`plot_fock_distribution`](@ref) | [Fock state](https://en.wikipedia.org/wiki/Fock_state) distribution |