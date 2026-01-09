#=
Functions for calculating metrics (distance measures) between states and operators.
=#

export fidelity
export tracedist, hilbert_dist, hellinger_dist
export bures_dist, bures_angle

@doc raw"""
    fidelity(ρ::QuantumObject, σ::QuantumObject)

Calculate the fidelity of two [`QuantumObject`](@ref):
``F(\hat{\rho}, \hat{\sigma}) = \textrm{Tr} \sqrt{\sqrt{\hat{\rho}} \hat{\sigma} \sqrt{\hat{\rho}}}``

Here, the definition is from [Nielsen-Chuang2011](@citet). It is the square root of the fidelity defined in [Jozsa1994](@citet).

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
function fidelity(ρ::QuantumObject{Operator}, σ::QuantumObject{Operator})
    sqrt_ρ = sqrt(ρ)
    eigval = abs.(eigvals(sqrt_ρ * σ * sqrt_ρ))
    return sum(sqrt, eigval)
end
fidelity(ρ::QuantumObject{Operator}, ψ::QuantumObject{Ket}) = sqrt(abs(expect(ρ, ψ)))
fidelity(ψ::QuantumObject{Ket}, σ::QuantumObject{Operator}) = fidelity(σ, ψ)
fidelity(ψ::QuantumObject{Ket}, ϕ::QuantumObject{Ket}) = abs(dot(ψ, ϕ))

@doc raw"""
    tracedist(ρ::QuantumObject, σ::QuantumObject)

Calculates the [trace distance](https://en.wikipedia.org/wiki/Trace_distance) between two [`QuantumObject`](@ref):
``T(\hat{\rho}, \hat{\sigma}) = \frac{1}{2} \lVert \hat{\rho} - \hat{\sigma} \rVert_1``

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
tracedist(
    ρ::QuantumObject{ObjType1},
    σ::QuantumObject{ObjType2},
) where {ObjType1 <: Union{Ket, Operator}, ObjType2 <: Union{Ket, Operator}} = norm(ket2dm(ρ) - ket2dm(σ), 1) / 2

@doc raw"""
    hilbert_dist(ρ::QuantumObject, σ::QuantumObject)

Calculates the Hilbert-Schmidt distance between two [`QuantumObject`](@ref):
``D_{HS}(\hat{\rho}, \hat{\sigma}) = \textrm{Tr}\left[\hat{A}^\dagger \hat{A}\right]``, where ``\hat{A} = \hat{\rho} - \hat{\sigma}``.

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).

# References
- [Vedral-Plenio1998](@citet)
"""
function hilbert_dist(
        ρ::QuantumObject{ObjType1},
        σ::QuantumObject{ObjType2},
    ) where {ObjType1 <: Union{Ket, Operator}, ObjType2 <: Union{Ket, Operator}}
    check_dimensions(ρ, σ)

    A = ket2dm(ρ) - ket2dm(σ)
    return tr(A' * A)
end

@doc raw"""
    hellinger_dist(ρ::QuantumObject, σ::QuantumObject)

Calculates the [Hellinger distance](https://en.wikipedia.org/wiki/Hellinger_distance) between two [`QuantumObject`](@ref):
``D_H(\hat{\rho}, \hat{\sigma}) = \sqrt{2 - 2 \textrm{Tr}\left(\sqrt{\hat{\rho}}\sqrt{\hat{\sigma}}\right)}``

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).

# References
- [Spehner2017](@citet)
"""
function hellinger_dist(
        ρ::QuantumObject{ObjType1},
        σ::QuantumObject{ObjType2},
    ) where {ObjType1 <: Union{Ket, Operator}, ObjType2 <: Union{Ket, Operator}}
    # Ket (pure state) doesn't need to do square root
    sqrt_ρ = isket(ρ) ? ket2dm(ρ) : sqrt(ρ)
    sqrt_σ = isket(σ) ? ket2dm(σ) : sqrt(σ)

    # `max` is to avoid numerical instabilities
    # it happens when ρ = σ, sum(eigvals) might be slightly larger than 1
    return sqrt(2.0 * max(0.0, 1.0 - sum(real, eigvals(sqrt_ρ * sqrt_σ))))
end

@doc raw"""
    bures_dist(ρ::QuantumObject, σ::QuantumObject)

Calculate the [Bures distance](https://en.wikipedia.org/wiki/Bures_metric) between two [`QuantumObject`](@ref):
``D_B(\hat{\rho}, \hat{\sigma}) = \sqrt{2 \left(1 - F(\hat{\rho}, \hat{\sigma}) \right)}``

Here, the definition of [`fidelity`](@ref) ``F`` is from [Nielsen-Chuang2011](@citet). It is the square root of the fidelity defined in [Jozsa1994](@citet).

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
bures_dist(ρ::QuantumObject, σ::QuantumObject) = sqrt(2 * (1 - fidelity(ρ, σ)))

@doc raw"""
    bures_angle(ρ::QuantumObject, σ::QuantumObject)

Calculate the [Bures angle](https://en.wikipedia.org/wiki/Bures_metric) between two [`QuantumObject`](@ref):
``D_A(\hat{\rho}, \hat{\sigma}) = \arccos\left(F(\hat{\rho}, \hat{\sigma})\right)``

Here, the definition of [`fidelity`](@ref) ``F`` is from [Nielsen-Chuang2011](@citet). It is the square root of the fidelity defined in [Jozsa1994](@citet).

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
bures_angle(ρ::QuantumObject, σ::QuantumObject) = acos(fidelity(ρ, σ))
