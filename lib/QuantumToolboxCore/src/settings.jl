Base.@kwdef mutable struct Settings
    tidyup_tol::Float64 = 1.0e-14
    auto_tidyup::Bool = true
    default_eltype_warning::Bool = true
end

#=
The settings below change a type, so they are stored as methods instead of fields: the compiler then knows their current
value, which keeps the code using them type-stable. Changing them redefines the method (see `Base.setproperty!` below).

Code that is already running keeps seeing the previous method (Julia's world age). To avoid silently creating quantum
objects with the previous element type, we store the world in which `default_eltype` was last changed, and throw an error
when it is read from an older world. A previous `progressmeter_kwargs` only changes how the progress bar looks, so it is
not checked.
=#
_default_eltype() = ComplexF64
_progressmeter_kwargs() = (showspeed = true, printed = true)

const _default_eltype_world = Ref(UInt(0))

# the world of the running code
# TODO: remove this branching and just use `Base.tls_world_age()` when Julia 1.10 is no longer supported
@static if isdefined(Base, :tls_world_age)
    _tls_world_age() = Base.tls_world_age()
else # Julia 1.10
    _tls_world_age() = ccall(:jl_get_tls_world_age, UInt, ())
end

@noinline _throw_previous_eltype_error() = throw(
    ErrorException(
        "`QuantumToolbox.settings.default_eltype` was changed by code that is still running, which would keep using the previous value. Change it with a separate statement, for example at the beginning of your script, and not inside a function or a block (such as `begin`, `let`, `try` or `for`) that also uses it.",
    ),
)

_depwarn_ProgressMeterKWARGS() = Base.depwarn(
    "`settings.ProgressMeterKWARGS` is deprecated and will be removed in next major release, use `settings.progressmeter_kwargs` instead.",
    :ProgressMeterKWARGS,
    force = true,
)

Base.@constprop :aggressive function Base.getproperty(s::Settings, name::Symbol)
    if name === :default_eltype
        (_tls_world_age() < _default_eltype_world[]) && _throw_previous_eltype_error()
        return _default_eltype()
    elseif name === :progressmeter_kwargs
        return _progressmeter_kwargs()
    elseif name === :ProgressMeterKWARGS
        _depwarn_ProgressMeterKWARGS()
        return getproperty(s, :progressmeter_kwargs)
    end
    return getfield(s, name)
end

function Base.setproperty!(s::Settings, name::Symbol, value)
    if name === :default_eltype
        (isconcretetype(value) && value <: FloatOrComplex) || throw(
            ArgumentError("`default_eltype` must be a concrete floating point or complex type, such as `ComplexF32` or `Float64`, but got `$value`."),
        )
        @eval _default_eltype() = $value
        _default_eltype_world[] = Base.get_world_counter()
    elseif name === :progressmeter_kwargs
        (value isa NamedTuple) || throw(ArgumentError("`progressmeter_kwargs` must be a `NamedTuple`, but got `$(typeof(value))`."))
        @eval _progressmeter_kwargs() = $value
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
QuantumToolbox.settings.default_eltype = ComplexF32
```

!!! warning "Change `progressmeter_kwargs` and `default_eltype` with a separate statement"
    Changing `progressmeter_kwargs` or `default_eltype` redefines a method, so code that is already running keeps seeing the previous value: reading `default_eltype` there throws an error, and the progress bar keeps the previous `progressmeter_kwargs`. Change them with a separate statement, for example at the beginning of your script, and not inside a function or a block (such as `begin`, `let`, `try` or `for`) that also uses them. For the same reason, they cannot be changed by a package while it is being precompiled: packages should pass the element type explicitly instead, e.g., `destroy(T, N)`.
"""
const settings = Settings()
