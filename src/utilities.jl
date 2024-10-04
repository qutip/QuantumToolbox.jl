#=
Utilities:
    internal (or external) functions which will be used throughout the entire package
=#

export gaussian, n_thermal, convert_unit
export row_major_reshape, meshgrid

@doc raw"""
    row_major_reshape(Q::AbstractArray, shapes...)

Reshapes `Q` in the row-major order, as numpy.
"""
row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} =
    PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

@doc raw"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    return X, Y
end

@doc raw"""
    gaussian(x::Number, μ::Number, σ::Number)

Returns the gaussian function ``\exp \left[- \frac{(x - \mu)^2}{2 \sigma^2} \right]``,
where ``\mu`` and ``\sigma^2`` are the mean and the variance respectively.
"""
gaussian(x::Number, μ::Number, σ::Number) = exp(-(x - μ)^2 / (2 * σ^2))

@doc raw"""
    n_thermal(ω::Real, ω_th::Real)

Return the number of photons in thermal equilibrium for an harmonic oscillator mode with frequency ``\omega``, at the temperature described by ``\omega_{\textrm{th}} \equiv k_B T / \hbar``:
```math
n(\omega, \omega_{\textrm{th}}) = \frac{1}{e^{\omega/\omega_{\textrm{th}}} - 1},
```
where ``\hbar`` is the reduced Planck constant, and ``k_B`` is the Boltzmann constant.
"""
function n_thermal(ω::T1, ω_th::T2) where {T1<:Real,T2<:Real}
    x = exp(ω / ω_th)
    n = ((x != 1) && (ω_th > 0)) ? 1 / (x - 1) : 0
    return _FType(promote_type(T1, T2))(n)
end

# some fundamental physical constants and common energy units
const _e = 1.602176565e-19 # elementary charge  (C)
const _kB = 1.3806488e-23  # Boltzmann constant (J/K)
const _h = 6.62607015e-34  # Planck constant    (J⋅s)
const _energy_units::Dict{Symbol,Float64} = Dict(
    # the values below are all in the unit of Joule
    :J => 1.0,
    :eV => _e,
    :meV => 1e-3 * _e,
    :GHz => 1e9 * _h,
    :mK => 1e-3 * _kB,
)

@doc raw"""
    convert_unit(value::Real, unit1::Symbol, unit2::Symbol)

Convert the energy `value` from `unit1` to `unit2`.

Note that `unit1` and `unit2` can be either the following `Symbol`:
- `:J`: Joule
- `:eV`: electron volt.
- `:meV`: milli-electron volt.
- `:GHz`: Giga-Hertz multiplied by Planck constant ``h``.
- `:mK`: milli-Kelvin multiplied by Boltzmann constant ``k_{\textrm{B}}``.

# Examples

```
julia> convert_unit(1, :eV, :J)
1.602176565e-19

julia> convert_unit(1, :GHz, :J)
6.62607015e-25

julia> convert_unit(1, :meV, :mK)
11604.51930280894
```
"""
function convert_unit(value::T, unit1::Symbol, unit2::Symbol) where {T<:Real}
    !haskey(_energy_units, unit1) && throw(ArgumentError("Invalid unit :$(unit1)"))
    !haskey(_energy_units, unit2) && throw(ArgumentError("Invalid unit :$(unit2)"))
    return _FType(T)(value * (_energy_units[unit1] / _energy_units[unit2]))
end

_get_dense_similar(A::AbstractArray, args...) = similar(A, args...)
_get_dense_similar(A::AbstractSparseMatrix, args...) = similar(nonzeros(A), args...)

_Ginibre_ensemble(n::Int, rank::Int = n) = randn(ComplexF64, n, rank) / sqrt(n)

makeVal(x::Val{T}) where {T} = x
makeVal(x) = Val(x)

getVal(x::Val{T}) where {T} = T

_get_size(A::AbstractMatrix) = size(A)
_get_size(A::AbstractVector) = (length(A), 1)

_non_static_array_warning(argname, arg::Tuple{}) =
    throw(ArgumentError("The argument $argname must be a Tuple or a StaticVector of non-zero length."))
_non_static_array_warning(argname, arg::Union{SVector{N,T},MVector{N,T},NTuple{N,T}}) where {N,T} = nothing
_non_static_array_warning(argname, arg::AbstractVector{T}) where {T} =
    @warn "The argument $argname should be a Tuple or a StaticVector for better performance. Try to use `$argname = $(Tuple(arg))` or `$argname = SVector(" *
          join(arg, ", ") *
          ")` instead of `$argname = $arg`." maxlog = 1

# functions for getting Float or Complex element type
_FType(::AbstractArray{T}) where {T<:Number} = _FType(T)
_FType(::Type{Int32}) = Float32
_FType(::Type{Int64}) = Float64
_FType(::Type{Float32}) = Float32
_FType(::Type{Float64}) = Float64
_FType(::Type{ComplexF32}) = Float32
_FType(::Type{ComplexF64}) = Float64
_CType(::AbstractArray{T}) where {T<:Number} = _CType(T)
_CType(::Type{Int32}) = ComplexF32
_CType(::Type{Int64}) = ComplexF64
_CType(::Type{Float32}) = ComplexF32
_CType(::Type{Float64}) = ComplexF64
_CType(::Type{ComplexF32}) = ComplexF32
_CType(::Type{ComplexF64}) = ComplexF64
