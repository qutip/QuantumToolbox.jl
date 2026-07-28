export PhysicalConstants, convert_unit

@doc raw"""
    const PhysicalConstants

A `NamedTuple` which stores some constant values listed in [*CODATA recommended values of the fundamental physical constants: 2022*](https://physics.nist.gov/cuu/pdf/wall_2022.pdf)

The current stored constants are:
- `c` : (exact) speed of light in vacuum with unit ``[\textrm{m}\cdot\textrm{s}^{-1}]``
- `G` : Newtonian constant of gravitation with unit ``[\textrm{m}^3\cdot\textrm{kg}^{−1}\cdot\textrm{s}^{−2}]``
- `h` : (exact) Planck constant with unit ``[\textrm{J}\cdot\textrm{s}]``
- `ħ` : reduced Planck constant with unit ``[\textrm{J}\cdot\textrm{s}]``
- `e` : (exact) elementary charge with unit ``[\textrm{C}]``
- `μ0` : vacuum magnetic permeability with unit ``[\textrm{N}\cdot\textrm{A}^{-2}]``
- `ϵ0` : vacuum electric permittivity with unit ``[\textrm{F}\cdot\textrm{m}^{-1}]``
- `k` : (exact) Boltzmann constant with unit ``[\textrm{J}\cdot\textrm{K}^{-1}]``
- `NA` : (exact) Avogadro constant with unit ``[\textrm{mol}^{-1}]``

# Examples

```jldoctest
julia> PhysicalConstants.ħ
1.0545718176461565e-34
```
"""
const PhysicalConstants = (
    c = 299792458.0,
    G = 6.6743e-11,
    h = 6.62607015e-34,
    ħ = 6.62607015e-34 / (2 * π),
    e = 1.602176634e-19,
    μ0 = 1.25663706127e-6,
    ϵ0 = 8.8541878188e-12,
    k = 1.380649e-23,
    NA = 6.02214076e23,
)

# common energy units (the values below are all in the unit of Joule)
const _energy_units::Dict{Symbol, Float64} = Dict(
    :J => 1.0,
    :eV => PhysicalConstants.e,
    :meV => 1.0e-3 * PhysicalConstants.e,
    :MHz => 1.0e6 * PhysicalConstants.h,
    :GHz => 1.0e9 * PhysicalConstants.h,
    :K => PhysicalConstants.k,
    :mK => 1.0e-3 * PhysicalConstants.k,
)

@doc raw"""
    convert_unit(value::Real, unit1::Symbol, unit2::Symbol)

Convert the energy `value` from `unit1` to `unit2`. The `unit1` and `unit2` can be either the following `Symbol`:
- `:J` : Joule
- `:eV` : electron volt
- `:meV` : milli-electron volt
- `:MHz` : Mega-Hertz multiplied by Planck constant ``h``
- `:GHz` : Giga-Hertz multiplied by Planck constant ``h``
- `:K` : Kelvin multiplied by Boltzmann constant ``k``
- `:mK` : milli-Kelvin multiplied by Boltzmann constant ``k``

Note that we use the values stored in [`PhysicalConstants`](@ref) to do the conversion.

# Examples

```jldoctest
julia> convert_unit(1, :eV, :J)
1.602176634e-19

julia> convert_unit(1, :GHz, :J)
6.62607015e-25

julia> round(convert_unit(1, :meV, :mK), digits=4)
11604.5181
```
"""
function convert_unit(value::T, unit1::Symbol, unit2::Symbol) where {T <: Real}
    !haskey(_energy_units, unit1) && throw(ArgumentError("Invalid unit :$(unit1)"))
    !haskey(_energy_units, unit2) && throw(ArgumentError("Invalid unit :$(unit2)"))
    return _float_type(T)(value * (_energy_units[unit1] / _energy_units[unit2]))
end
