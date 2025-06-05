Base.@kwdef mutable struct Settings
    tidyup_tol::Float64 = 1e-14
    auto_tidyup::Bool = true
end

function Base.show(io::IO, s::Settings)
    # To align the output and make it easier to read
    # we use rpad `11`, which is the length of string: `auto_tidyup`
    println(io, "QuantumToolbox.jl Settings")
    println(io, "--------------------------")
    map(n -> println(io, rpad("$n", 11, " "), " = ", getfield(s, n)), fieldnames(Settings))
    return nothing
end

@doc raw"""
    QuantumToolbox.settings 

Contains all the default global settings of QuantumToolbox.jl.

# List of settings

- `tidyup_tol::Float64 = 1e-14` : tolerance for [`tidyup`](@ref) and [`tidyup!`](@ref).
- `auto_tidyup::Bool = true` : Automatically tidyup.

For detailed explanation of each settings, see our documentation [here](https://qutip.org/QuantumToolbox.jl/stable/users_guide/settings).

# Change default settings

One can overwrite the default global settings by

```julia
using QuantumToolbox

QuantumToolbox.settings.tidyup_tol = 1e-10
QuantumToolbox.settings.auto_tidyup = false
```
"""
const settings = Settings()
