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
- `progressmeter_kwargs::NamedTuple = (showspeed = true, printed = true)` : Default keyword arguments for progress bar in [`ProgressMeter.jl`](https://github.com/timholy/ProgressMeter.jl). This allows the customization of progress bar.
- `default_eltype::Type = ComplexF64` : Element type of the quantum objects created without specifying it, such as `destroy(N)` or `fock(N, j)`. It must be a concrete floating point or complex type, such as `ComplexF32` or `Float64`. A real type can be used for real-valued systems: functions that need complex numbers, such as [`sigmay`](@ref) or [`rand_ket`](@ref), then use its complex version.
- `show_eltype_warning::Bool = true` : Warn (once) when a quantum object is created with an element type more general than `default_eltype`, e.g., `0.5 * destroy(N)` with `default_eltype = ComplexF32`, or a complex state with a real `default_eltype`.

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

## Change the default element type

By default, quantum objects are created with element type `ComplexF64`. To use a lower precision, for example to run faster or to save memory, one can change `default_eltype`:

```@example settings
QuantumToolbox.settings.default_eltype = ComplexF32

a = destroy(5)
eltype(a)
```

Note that a `Float64` number promotes the element type back to `ComplexF64`, following the usual promotion rules of Julia. A warning is shown (once per session) when this happens:

```@example settings
eltype(0.5 * a)
```

To keep the lower precision, use `Float32` numbers instead, such as `0.5f0`:

```@example settings
eltype(0.5f0 * a)
```

An element type given explicitly is still used. The warning above cannot tell it apart from a promotion, so set `QuantumToolbox.settings.show_eltype_warning = false` if a more general type is intended:

```@example settings
eltype(destroy(ComplexF64, 5))
```

```@setup settings
# restore the default settings, since the other pages of the documentation run in the same Julia session
QuantumToolbox.settings.tidyup_tol = 1.0e-14
QuantumToolbox.settings.auto_tidyup = true
QuantumToolbox.settings.default_eltype = ComplexF64
```

!!! warning "Change `progressmeter_kwargs` and `default_eltype` with a separate statement"
    Changing `progressmeter_kwargs` or `default_eltype` redefines a method, so code that is already running keeps seeing the previous value: reading `default_eltype` there throws an error, and the progress bar keeps the previous `progressmeter_kwargs`. Change them with a separate statement, for example at the beginning of your script, and not inside a function or a block (such as `begin`, `let`, `try` or `for`) that also uses them. For the same reason, they cannot be changed by a package while it is being precompiled: packages should pass the element type explicitly instead, e.g., `destroy(T, N)`.
