#=
Functions for calculating metrics (distance measures) between states and operators.
=#

export entropy_vn, entanglement, tracedist, fidelity

@doc raw"""
    entropy_vn(ρ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy)
``S = - \Tr \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right]`` where ``\hat{\rho}``
is the density matrix of the system.

The `base` parameter specifies the base of the logarithm to use, and when using the default value 0,
the natural logarithm is used. The `tol` parameter
describes the absolute tolerance for detecting the zero-valued eigenvalues of the density
matrix ``\hat{\rho}``.

# Examples

Pure state:
```
julia> ψ = fock(2,0)
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im

julia> ρ = ket2dm(ψ)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> entropy_vn(ρ, base=2)
-0.0
```

Mixed state:
```
julia> ρ = 1 / 2 * ( ket2dm(fock(2,0)) + ket2dm(fock(2,1)) )
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im

julia> entropy_vn(ρ, base=2)
1.0
```
"""
function entropy_vn(
    ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject};
    base::Int = 0,
    tol::Real = 1e-15,
) where {T}
    vals = eigvals(ρ)
    indexes = abs.(vals) .> tol
    1 ∉ indexes && return 0
    nzvals = vals[indexes]
    logvals = base != 0 ? log.(base, Complex.(nzvals)) : log.(Complex.(nzvals))
    return -real(sum(nzvals .* logvals))
end

"""
    entanglement(QO::QuantumObject, sel::Vector)

Calculates the entanglement by doing the partial trace of `QO`, selecting only the dimensions
with the indices contained in the `sel` vector, and then using the Von Neumann entropy [`entropy_vn`](@ref).
"""
function entanglement(
    QO::QuantumObject{<:AbstractArray{T},OpType},
    sel::Vector{Int},
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    ψ = normalize(QO)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end
entanglement(QO::QuantumObject, sel::Int) = entanglement(QO, [sel])

@doc raw"""
    tracedist(ρ::QuantumObject, σ::QuantumObject)

Calculates the [trace distance](https://en.wikipedia.org/wiki/Trace_distance) between two [`QuantumObject`](@ref):
``T(\rho, \sigma) = \frac{1}{2} \lVert \rho - \sigma \rVert_1``

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
tracedist(
    ρ::QuantumObject{<:AbstractArray{T1},ObjType1},
    σ::QuantumObject{<:AbstractArray{T2},ObjType2},
) where {
    T1,
    T2,
    ObjType1<:Union{KetQuantumObject,OperatorQuantumObject},
    ObjType2<:Union{KetQuantumObject,OperatorQuantumObject},
} = norm(ket2dm(ρ) - ket2dm(σ), 1) / 2

@doc raw"""
    fidelity(ρ::QuantumObject, σ::QuantumObject)

Calculate the fidelity of two [`QuantumObject`](@ref):
``F(\rho, \sigma) = \textrm{Tr} \sqrt{\sqrt{\rho} \sigma \sqrt{\rho}}``

Here, the definition is from Nielsen & Chuang, "Quantum Computation and Quantum Information". It is the square root of the fidelity defined in R. Jozsa, Journal of Modern Optics, 41:12, 2315 (1994).

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
function fidelity(
    ρ::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    σ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    sqrt_ρ = sqrt(ρ)
    eigval = abs.(eigvals(sqrt_ρ * σ * sqrt_ρ))
    return sum(sqrt, eigval)
end
fidelity(
    ρ::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2} = sqrt(abs(expect(ρ, ψ)))
fidelity(
    ψ::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    σ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2} = fidelity(σ, ψ)
fidelity(
    ψ::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    ϕ::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2} = abs(dot(ψ, ϕ))
