Base.@kwdef mutable struct Settings
    tidyup_tol::Float64 = 1.0e-14
    auto_tidyup::Bool = true
    default_eltype_warning::Bool = true
end

#=
The settings below change a type, so they are stored as methods instead of fields: the compiler then knows their current
value, which keeps the code using them type-stable. Changing them redefines the method (see `Base.setproperty!` below).

Code that is already running keeps seeing the previous method (Julia's world age). To avoid silently using the previous
value, we store the world in which each setting was last changed, and throw an error when it is read from an older world.
=#
_default_eltype() = ComplexF64
_progressmeter_kwargs() = (showspeed = true, printed = true)

const _default_eltype_world = Ref(UInt(0))
const _progressmeter_kwargs_world = Ref(UInt(0))

function _check_world(changed_world::Base.RefValue{UInt}, name::Symbol)
    (ccall(:jl_get_tls_world_age, UInt, ()) < changed_world[]) && throw(
        ErrorException(
            "`QuantumToolbox.settings.$name` was changed inside code that is still running, which would keep using the previous value. Change it outside of that code, for example at the beginning of your script.",
        ),
    )
    return nothing
end

# A package being precompiled cannot redefine methods of another package, because the change could not be saved
function _check_not_precompiling(name::Symbol)
    (ccall(:jl_generating_output, Cint, ()) == 1) && throw(
        ErrorException(
            "`QuantumToolbox.settings.$name` cannot be changed while a package is being precompiled. Packages should not change it for their users: pass the element type explicitly instead, e.g., `destroy(T, N)`.",
        ),
    )
    return nothing
end

_depwarn_ProgressMeterKWARGS() = Base.depwarn(
    "`settings.ProgressMeterKWARGS` is deprecated and will be removed in next major release, use `settings.progressmeter_kwargs` instead.",
    :ProgressMeterKWARGS,
    force = true,
)

Base.@constprop :aggressive function Base.getproperty(s::Settings, name::Symbol)
    if name === :default_eltype
        _check_world(_default_eltype_world, name)
        return _default_eltype()
    elseif name === :progressmeter_kwargs
        _check_world(_progressmeter_kwargs_world, name)
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
        _check_not_precompiling(name)
        @eval _default_eltype() = $value
        _default_eltype_world[] = Base.get_world_counter()
    elseif name === :progressmeter_kwargs
        (value isa NamedTuple) || throw(ArgumentError("`progressmeter_kwargs` must be a `NamedTuple`, but got `$(typeof(value))`."))
        _check_not_precompiling(name)
        @eval _progressmeter_kwargs() = $value
        _progressmeter_kwargs_world[] = Base.get_world_counter()
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

#=
Warn (once) when a quantum object has a higher precision than `default_eltype`. This only applies when the default has a
lower precision than `Float64`, because only then a `Float64` number can silently promote it (e.g., `1.0 * destroy(N)` with
`default_eltype = ComplexF32`). With higher defaults, a higher precision can only be chosen explicitly.
=#
function _has_higher_precision(::Type{T}) where {T <: AbstractFloat}
    D = real(_default_eltype())
    return (promote_type(D, Float64) !== D) && (promote_type(T, D) !== D)
end
_has_higher_precision(::Type{Complex{T}}) where {T <: AbstractFloat} = _has_higher_precision(T)
_has_higher_precision(::Type) = false

function _check_eltype(::Type{T}) where {T}
    if _has_higher_precision(T) && settings.default_eltype_warning
        @warn "A quantum object with element type `$T` was created, which has a higher precision than `QuantumToolbox.settings.default_eltype = $(_default_eltype())`. This usually happens when a `Float64` number, such as `0.5` or `2π`, multiplies a quantum object. This warning is shown only once. To disable it, set `QuantumToolbox.settings.default_eltype_warning = false`." maxlog =
            1
    end
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
- `default_eltype_warning::Bool = true` : Warn (once) when a quantum object is created with a higher precision than `default_eltype`. This only applies when `default_eltype` has a lower precision than `Float64`.

For detailed explanation of each settings, see our documentation [here](https://qutip.org/QuantumToolbox.jl/stable/users_guide/settings).

# Change default settings

One can overwrite the default global settings by

```julia
using QuantumToolbox

QuantumToolbox.settings.tidyup_tol = 1e-10
QuantumToolbox.settings.auto_tidyup = false
QuantumToolbox.settings.default_eltype = ComplexF32
```

!!! warning "Change `progressmeter_kwargs` and `default_eltype` outside of the code using them"
    Changing `progressmeter_kwargs` or `default_eltype` redefines a method, so code that is already running keeps using the previous value. If this code reads the setting, an error is thrown. Change them outside of that code, for example at the beginning of your script. For the same reason, they cannot be changed by a package while it is being precompiled: packages should pass the element type explicitly instead, e.g., `destroy(T, N)`.
"""
const settings = Settings()
