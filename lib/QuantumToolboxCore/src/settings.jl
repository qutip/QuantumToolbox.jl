Base.@kwdef mutable struct Settings
    tidyup_tol::Float64 = 1.0e-14
    auto_tidyup::Bool = true
    default_eltype_warning::Bool = true
end

#=
The settings below change a type, so they are stored as compile-time preferences (see Preferences.jl): their values are
constants, which keeps the code using them type-stable. They are saved as Julia code (e.g., "ComplexF32") in the
`LocalPreferences.toml` file of the active project, and are read when QuantumToolboxCore is precompiled. So, a change is
used after restarting Julia, which then precompiles QuantumToolboxCore and the packages depending on it again.
=#
function _load_preference(name::String, default::String, isvalid::Function)
    code = load_preference(@__MODULE__, name, default)
    value = try
        eval(Meta.parse(code))
    catch
        nothing
    end
    isvalid(value) ||
        error("Invalid preference `$name = $code` for QuantumToolboxCore. Fix or remove it in the `LocalPreferences.toml` file of your project.")
    return value
end

function _save_preference(name::String, value)
    code = repr(value)
    value_read = try
        eval(Meta.parse(code))
    catch
        nothing
    end
    isequal(value_read, value) || throw(
        ArgumentError(
            "`$name = $code` cannot be saved as a preference, because it cannot be read back when QuantumToolboxCore is precompiled. Only values written with the names defined in Julia's `Base` can be saved, e.g., `ComplexF32`, `Complex{BigFloat}` or `(showspeed = false, output = devnull)`.",
        ),
    )
    set_preferences!(@__MODULE__, name => code; force = true)
    @info "`QuantumToolbox.settings.$name = $code` was saved in the `LocalPreferences.toml` file of your project. Restart Julia to use it."
    return nothing
end

_isvalid_eltype(T) = T isa Type && isconcretetype(T) && T <: FloatOrComplex

const _DEFAULT_ELTYPE = _load_preference("default_eltype", "ComplexF64", _isvalid_eltype)
const _PROGRESSMETER_KWARGS = _load_preference("progressmeter_kwargs", "(showspeed = true, printed = true)", kw -> kw isa NamedTuple)

_depwarn_ProgressMeterKWARGS() = Base.depwarn(
    "`settings.ProgressMeterKWARGS` is deprecated and will be removed in next major release, use `settings.progressmeter_kwargs` instead.",
    :ProgressMeterKWARGS,
    force = true,
)

Base.@constprop :aggressive function Base.getproperty(s::Settings, name::Symbol)
    if name === :default_eltype
        return _DEFAULT_ELTYPE
    elseif name === :progressmeter_kwargs
        return _PROGRESSMETER_KWARGS
    elseif name === :ProgressMeterKWARGS
        _depwarn_ProgressMeterKWARGS()
        return _PROGRESSMETER_KWARGS
    end
    return getfield(s, name)
end

function Base.setproperty!(s::Settings, name::Symbol, value)
    if name === :default_eltype
        _isvalid_eltype(value) || throw(
            ArgumentError("`default_eltype` must be a concrete floating point or complex type, such as `ComplexF32` or `Float64`, but got `$value`."),
        )
        _save_preference("default_eltype", value)
    elseif name === :progressmeter_kwargs
        (value isa NamedTuple) || throw(ArgumentError("`progressmeter_kwargs` must be a `NamedTuple`, but got `$(typeof(value))`."))
        _save_preference("progressmeter_kwargs", value)
    elseif name === :ProgressMeterKWARGS
        _depwarn_ProgressMeterKWARGS()
        setproperty!(s, :progressmeter_kwargs, value)
    else
        setfield!(s, name, convert(fieldtype(Settings, name), value))
    end
    return value
end

Base.propertynames(::Settings, ::Bool = false) =
    (:tidyup_tol, :auto_tidyup, :progressmeter_kwargs, :default_eltype, :default_eltype_warning)

function Base.show(io::IO, s::Settings)
    # get maximum property name to align the output and make it easier to read
    maxLen = maximum(length ∘ String, propertynames(s))

    println(io, "QuantumToolbox.jl Settings")
    println(io, "--------------------------")
    map(n -> println(io, rpad("$n", maxLen, " "), " = ", getproperty(s, n)), propertynames(s))
    return nothing
end

@doc raw"""
    QuantumToolbox.settings

Contains all the default global settings of QuantumToolbox.jl and all related libraries.

# List of settings

- `tidyup_tol::Float64 = 1.0e-14` : tolerance for [`tidyup`](@ref QuantumToolboxCore.tidyup), [`tidyup!`](@ref QuantumToolboxCore.tidyup!) and [`to_sparse`](@ref QuantumToolboxCore.to_sparse).
- `auto_tidyup::Bool = true` : Automatically `tidyup`.
- `progressmeter_kwargs::NamedTuple = (showspeed = true, printed = true)` : Default keyword arguments for progress bar in [`ProgressMeter.jl`](https://github.com/timholy/ProgressMeter.jl). This allows the customization of progress bar.
- `default_eltype::Type = ComplexF64` : Element type of the quantum objects created without specifying it, e.g., `destroy(N)` or `fock(N, j)`. It must be a concrete floating point or complex type, such as `ComplexF32` or `Float64`.
- `default_eltype_warning::Bool = true` : Warn (once) when a quantum object is created with an element type more general than `default_eltype`, e.g., `0.5 * destroy(N)` with `default_eltype = ComplexF32`.

For detailed explanation of each settings, see our documentation [here](https://qutip.org/QuantumToolbox.jl/stable/users_guide/settings).

# Change default settings

One can overwrite the default global settings by

```julia
using QuantumToolbox

QuantumToolbox.settings.tidyup_tol = 1e-10
QuantumToolbox.settings.auto_tidyup = false
QuantumToolbox.settings.default_eltype = ComplexF32 # used after restarting Julia
```

!!! note "`progressmeter_kwargs` and `default_eltype` are used after restarting Julia"
    These two settings change a type, so they are saved as preferences (using [`Preferences.jl`](https://github.com/JuliaPackaging/Preferences.jl)) in the `LocalPreferences.toml` file of the active project, and are used after restarting Julia. Only values written with the names defined in Julia's `Base` can be saved, e.g., `ComplexF32`, `Complex{BigFloat}` or `(showspeed = false, output = devnull)`.
"""
const settings = Settings()
