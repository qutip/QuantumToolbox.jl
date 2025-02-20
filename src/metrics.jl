#=
Functions for calculating metrics (distance measures) between states and operators.
=#

export fidelity
export tracedist, hilbert_dist
export bures_dist, bures_angle

@doc raw"""
    fidelity(ρ::QuantumObject, σ::QuantumObject)

Calculate the fidelity of two [`QuantumObject`](@ref):
``F(\hat{\rho}, \hat{\sigma}) = \textrm{Tr} \sqrt{\sqrt{\hat{\rho}} \hat{\sigma} \sqrt{\hat{\rho}}}``

Here, the definition is from [Nielsen-Chuang2011](@citet). It is the square root of the fidelity defined in [Jozsa1994](@citet).

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
function fidelity(ρ::QuantumObject{OperatorQuantumObject}, σ::QuantumObject{OperatorQuantumObject})
    sqrt_ρ = sqrt(ρ)
    eigval = abs.(eigvals(sqrt_ρ * σ * sqrt_ρ))
    return sum(sqrt, eigval)
end
fidelity(ρ::QuantumObject{OperatorQuantumObject}, ψ::QuantumObject{KetQuantumObject}) = sqrt(abs(expect(ρ, ψ)))
fidelity(ψ::QuantumObject{KetQuantumObject}, σ::QuantumObject{OperatorQuantumObject}) = fidelity(σ, ψ)
fidelity(ψ::QuantumObject{KetQuantumObject}, ϕ::QuantumObject{KetQuantumObject}) = abs(dot(ψ, ϕ))

@doc raw"""
    tracedist(ρ::QuantumObject, σ::QuantumObject)

Calculates the [trace distance](https://en.wikipedia.org/wiki/Trace_distance) between two [`QuantumObject`](@ref):
``T(\hat{\rho}, \hat{\sigma}) = \frac{1}{2} \lVert \hat{\rho} - \hat{\sigma} \rVert_1``

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
tracedist(
    ρ::QuantumObject{ObjType1},
    σ::QuantumObject{ObjType2},
) where {
    ObjType1<:Union{KetQuantumObject,OperatorQuantumObject},
    ObjType2<:Union{KetQuantumObject,OperatorQuantumObject},
} = norm(ket2dm(ρ) - ket2dm(σ), 1) / 2

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
) where {
    ObjType1<:Union{KetQuantumObject,OperatorQuantumObject},
    ObjType2<:Union{KetQuantumObject,OperatorQuantumObject},
}
    check_dimensions(ρ, σ)

    A = ket2dm(ρ) - ket2dm(σ)
    return tr(A' * A)
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
