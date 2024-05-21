# [Installation](@id doc: Installation)

## QuantumToolbox.jl
To install `QuantumToolbox.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("QuantumToolbox")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.7) pkg> add QuantumToolbox
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).  
!!! note "Julia 1.7"
    `QuantumToolbox.jl` requires Julia 1.7 or higher.

## Verify the installation
To load the package and check the version information, use either `versioninfo()` or `about()`, namely
```julia
julia> using QuantumToolbox
julia> QuantumToolbox.versioninfo()
julia> QuantumToolbox.about()
```