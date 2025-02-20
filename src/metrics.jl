#=
Functions for calculating metrics (distance measures) between states and operators.
=#

export fidelity
export tracedist, hilbert_dist

@doc raw"""
    fidelity(ρ::QuantumObject, σ::QuantumObject)

Calculate the fidelity of two [`QuantumObject`](@ref):
``F(\hat{\rho}, \hat{\sigma}) = \textrm{Tr} \sqrt{\sqrt{\hat{\rho}} \hat{\sigma} \sqrt{\hat{\rho}}}``

Here, the definition is from Nielsen & Chuang, "Quantum Computation and Quantum Information". It is the square root of the fidelity defined in R. Jozsa, Journal of Modern Optics, 41:12, 2315 (1994).

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
