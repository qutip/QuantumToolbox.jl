#=
Functions for calculating metrics (distance measures) between states and operators.
=#

export tracedist, fidelity

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
