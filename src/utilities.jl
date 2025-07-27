#=
Utilities:
    internal (or external) functions which will be used throughout the entire package
=#

export gaussian, n_thermal
export PhysicalConstants, convert_unit
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
    return _float_type(promote_type(T1, T2))(n)
end

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
    G = 6.67430e-11,
    h = 6.62607015e-34,
    ħ = 6.62607015e-34 / (2 * π),
    e = 1.602176634e-19,
    μ0 = 1.25663706127e-6,
    ϵ0 = 8.8541878188e-12,
    k = 1.380649e-23,
    NA = 6.02214076e23,
)

# common energy units (the values below are all in the unit of Joule)
const _energy_units::Dict{Symbol,Float64} = Dict(
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
function convert_unit(value::T, unit1::Symbol, unit2::Symbol) where {T<:Real}
    !haskey(_energy_units, unit1) && throw(ArgumentError("Invalid unit :$(unit1)"))
    !haskey(_energy_units, unit2) && throw(ArgumentError("Invalid unit :$(unit2)"))
    return _float_type(T)(value * (_energy_units[unit1] / _energy_units[unit2]))
end

get_typename_wrapper(A) = Base.typename(typeof(A)).wrapper

_dense_similar(A::AbstractArray, args...) = similar(A, args...)
_dense_similar(A::AbstractSparseMatrix, args...) = similar(nonzeros(A), args...)

_sparse_similar(A::AbstractArray, args...) = sparse(args...)

_Ginibre_ensemble(n::Int, rank::Int = n) = randn(ComplexF64, n, rank) / sqrt(n)

makeVal(x::Val{T}) where {T} = x
makeVal(x) = Val(x)

getVal(x::Val{T}) where {T} = T
getVal(x) = x # getVal for any other type

_get_size(A::AbstractMatrix) = size(A)
_get_size(A::AbstractVector) = (length(A), 1)
_get_size(A::AbstractSciMLOperator) = size(A)

_non_static_array_warning(argname, arg::Tuple{}) =
    throw(ArgumentError("The argument $argname must be a Tuple or a StaticVector of non-zero length."))
_non_static_array_warning(argname, arg::Union{SVector{N,T},MVector{N,T},NTuple{N,T}}) where {N,T} = nothing
_non_static_array_warning(argname, arg::AbstractVector{T}) where {T} =
    @warn "The argument $argname should be a Tuple or a StaticVector for better performance. Try to use `$argname = $(Tuple(arg))` instead of `$argname = $arg`. " *
          "Alternatively, you can do `import QuantumToolbox: SVector` " *
          "and use `$argname = SVector(" *
          join(arg, ", ") *
          ")`." maxlog = 1

# lazy tensor warning
for AType in (:AbstractArray, :AbstractSciMLOperator)
    for BType in (:AbstractArray, :AbstractSciMLOperator)
        if AType == BType == :AbstractArray
            @eval begin
                _lazy_tensor_warning(::$AType, ::$BType) = nothing
            end
        else
            @eval begin
                _lazy_tensor_warning(A::$AType, B::$BType) =
                    @warn "using lazy tensor (which can hurt performance) between data types: $(get_typename_wrapper(A)) and $(get_typename_wrapper(B))"
            end
        end
    end
end

# functions for getting Float or Complex element type
_float_type(::AbstractArray{T}) where {T<:Number} = _float_type(T)
_float_type(::Type{Int32}) = Float32
_float_type(::Type{Int64}) = Float64
_float_type(::Type{Float32}) = Float32
_float_type(::Type{Float64}) = Float64
_float_type(::Type{Complex{Int32}}) = Float32
_float_type(::Type{Complex{Int64}}) = Float64
_float_type(::Type{Complex{Float32}}) = Float32
_float_type(::Type{Complex{Float64}}) = Float64
_float_type(T::Type{<:Real}) = T # Allow other untracked Real types, like ForwardDiff.Dual
_complex_float_type(::AbstractArray{T}) where {T<:Number} = _complex_float_type(T)
_complex_float_type(::Type{Int32}) = ComplexF32
_complex_float_type(::Type{Int64}) = ComplexF64
_complex_float_type(::Type{Float32}) = ComplexF32
_complex_float_type(::Type{Float64}) = ComplexF64
_complex_float_type(::Type{Complex{Int32}}) = ComplexF32
_complex_float_type(::Type{Complex{Int64}}) = ComplexF64
_complex_float_type(::Type{Complex{Float32}}) = ComplexF32
_complex_float_type(::Type{Complex{Float64}}) = ComplexF64
_complex_float_type(T::Type{<:Complex}) = T # Allow other untracked Complex types, like ForwardDiff.Dual

_convert_eltype_wordsize(::Type{T}, ::Val{64}) where {T<:Int} = Int64
_convert_eltype_wordsize(::Type{T}, ::Val{32}) where {T<:Int} = Int32
_convert_eltype_wordsize(::Type{T}, ::Val{64}) where {T<:AbstractFloat} = Float64
_convert_eltype_wordsize(::Type{T}, ::Val{32}) where {T<:AbstractFloat} = Float32
_convert_eltype_wordsize(::Type{Complex{T}}, ::Val{64}) where {T<:Union{Int,AbstractFloat}} = ComplexF64
_convert_eltype_wordsize(::Type{Complex{T}}, ::Val{32}) where {T<:Union{Int,AbstractFloat}} = ComplexF32
