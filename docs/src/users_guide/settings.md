# [QuantumToolbox Settings](@id doc:QuantumToolbox-Settings)

In this section, we introduce the default global settings used throughout the package and show how to modify them.

All settings are stored in [`QuantumToolbox.settings`](@ref).

!!! warning "Differences from QuTiP"
    Due to the differences in programming languages, solving algorithms, and many other reasons, these global settings (including their default values and usage) may be very different from those in `Python QuTiP`.

## List of settings

Here, we list out each setting along with the specific functions that will use it.

- `tidyup_tol::Real = 1e-14` : tolerance for [`tidyup`](@ref) and [`tidyup!`](@ref).
- `auto_tidyup::Bool = true` : Automatically tidyup during the following situations:
    * Solving for eigenstates, including [`eigenstates`](@ref), [`eigsolve`](@ref), and [`eigsolve_al`](@ref).
- (to be announced)

## Change default settings

First, we can check the current [`QuantumToolbox.settings`](@ref):

```@example settings
using QuantumToolbox

QuantumToolbox.settings
```

Next, one can overwrite the default settings by

```@example settings
QuantumToolbox.settings.tidyup_tol[] = 1e-10
QuantumToolbox.settings.auto_tidyup[] = false

QuantumToolbox.settings
```