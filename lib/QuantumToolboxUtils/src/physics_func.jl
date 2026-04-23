export gaussian, n_thermal

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
function n_thermal(ω::T1, ω_th::T2) where {T1 <: Real, T2 <: Real}
    x = exp(ω / ω_th)
    n = ((x != 1) && (ω_th > 0)) ? 1 / (x - 1) : 0
    return _float_type(promote_type(T1, T2))(n)
end

_Ginibre_ensemble(rng::AbstractRNG, ::Type{T}, n::Int, rank::Int = n) where {T <: Complex} = randn(rng, T, n, rank) / sqrt(T(n))

_Boltzmann_weight(β::T, E::Int) where {T <: Real} = (E != 0 || isfinite(β)) ? exp(-β * E) : one(T)
