# [QuantumToolbox Settings](@id doc:QuantumToolbox-Settings)

In this section, we introduce the default global settings used throughout the package and show how to modify them.

All settings are stored in [`QuantumToolbox.settings`](@ref QuantumToolboxCore.settings).

!!! warning "Differences from QuTiP"
    Due to the differences in programming languages, solving algorithms, and many other reasons, these global settings (including their default values and usage) may be very different from those in `Python QuTiP`.

## List of settings

Here, we list out each setting along with the specific functions that will use it.

- `tidyup_tol::Float64 = 1.0e-14` : tolerance for [`tidyup`](@ref), [`tidyup!`](@ref) and [`to_sparse`](@ref).
- `auto_tidyup::Bool = true` : Automatically tidyup during the following situations:
    * Solving for eigenstates, including [`eigenstates`](@ref), [`eigsolve`](@ref), [`eigsolve_al`](@ref)
    * Creating [`bloch_redfield_tensor`](@ref) or [`brterm`](@ref), and solving [`brmesolve`](@ref).
- `progressmeter_kwargs::NamedTuple = (showspeed = true, printed = true)` : Default keyword arguments for progress bar in [`ProgressMeter.jl`](https://github.com/timholy/ProgressMeter.jl). This allows the customization of progress bar. It is saved as a preference, see [below](@ref doc:QuantumToolbox-Settings:preferences).
- `default_eltype::Type = ComplexF64` : Element type of the quantum objects created without specifying it, such as `destroy(N)` or `fock(N, j)`. It must be a concrete floating point or complex type, such as `ComplexF32` or `Float64`. A real type can be used for real-valued systems: functions that need complex numbers, such as [`sigmay`](@ref) or [`rand_ket`](@ref), then use its complex version. It is saved as a preference, see [below](@ref doc:QuantumToolbox-Settings:preferences).
- `default_eltype_warning::Bool = true` : Warn (once) when a quantum object is created with an element type more general than `default_eltype`, e.g., `0.5 * destroy(N)` with `default_eltype = ComplexF32`, or a complex state with a real `default_eltype`.

## Change default settings

First, we can check the current [`QuantumToolbox.settings`](@ref QuantumToolboxCore.settings):

```@example settings
using QuantumToolbox

QuantumToolbox.settings
```

Next, one can overwrite the default settings by

```@example settings
QuantumToolbox.settings.tidyup_tol = 1e-10
QuantumToolbox.settings.auto_tidyup = false

QuantumToolbox.settings
```

```@setup settings
# restore the default settings, since the other pages of the documentation run in the same Julia session
QuantumToolbox.settings.tidyup_tol = 1.0e-14
QuantumToolbox.settings.auto_tidyup = true
```

## [Settings saved as preferences](@id doc:QuantumToolbox-Settings:preferences)

`default_eltype` and `progressmeter_kwargs` change the type of the objects created by `QuantumToolbox.jl`. To keep the code using them type-stable, they are saved as preferences (using [`Preferences.jl`](https://github.com/JuliaPackaging/Preferences.jl)), which are read when `QuantumToolbox.jl` is precompiled. So, a change is used only after restarting Julia:

```julia
using QuantumToolbox

QuantumToolbox.settings.default_eltype = ComplexF32
QuantumToolbox.settings.progressmeter_kwargs = (showspeed = false, output = devnull)
```

This saves them in the `LocalPreferences.toml` file of the active project (next to its `Project.toml`), which can also be edited directly:

```toml
[QuantumToolboxCore]
default_eltype = "ComplexF32"
progressmeter_kwargs = "(showspeed = false, output = Base.DevNull())"
```

After restarting Julia, `QuantumToolbox.jl` and the packages depending on it are precompiled again, and the new values are used:

```julia
using QuantumToolbox

a = destroy(5)
eltype(a)         # ComplexF32
```

Note that a `Float64` number promotes the element type back to `ComplexF64`, following the usual promotion rules of Julia, and a warning is shown (once per session). To keep the lower precision, use `Float32` numbers instead, such as `0.5f0`:

```julia
eltype(0.5 * a)   # ComplexF64, with a warning
eltype(0.5f0 * a) # ComplexF32
```

An element type given explicitly, such as `destroy(ComplexF64, 5)`, is still used. The warning cannot tell it apart from a promotion, so set `QuantumToolbox.settings.default_eltype_warning = false` if a more general type is intended.

!!! note "Which values can be saved?"
    The preferences are saved as Julia code, which is read back when `QuantumToolbox.jl` is precompiled. So, only values written with the names defined in Julia's `Base` can be saved, e.g., `ComplexF32`, `Complex{BigFloat}` or `(showspeed = false, output = devnull)`, but not `Complex{Double64}` (from `DoubleFloats.jl`) or `output = stdout`. To go back to the default values, set them again (e.g., `QuantumToolbox.settings.default_eltype = ComplexF64`) or remove them from `LocalPreferences.toml`. If the active project does not depend directly on `QuantumToolboxCore` (the library defining the settings), it is also added to the `[extras]` of its `Project.toml`.
