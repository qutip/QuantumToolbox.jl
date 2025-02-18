#=
Entropy related functions.
=#

export entropy_vn
export entanglement

@doc raw"""
    entropy_vn(ρ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy) ``S = - \textrm{Tr} \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right]``, where ``\hat{\rho}`` is the density matrix of the system.

# Notes

- `ρ` is the quantum state, can be either a [`Ket`](@ref) or [`Operator`](@ref).
- `base` specifies the base of the logarithm to use, and when using the default value `0`, the natural logarithm is used.
- `tol` describes the absolute tolerance for detecting the zero-valued eigenvalues of the density matrix ``\hat{\rho}``.

# Examples

Pure state:
```jldoctest
julia> ψ = fock(2,0)

Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im

julia> entropy_vn(ψ, base=2)
-0.0
```

Mixed state:
```jldoctest
julia> ρ = maximally_mixed_dm(2)

Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Diagonal{ComplexF64, Vector{ComplexF64}}:
 0.5-0.0im      ⋅    
     ⋅      0.5-0.0im

julia> entropy_vn(ρ, base=2)
1.0
```
"""
function entropy_vn(
    ρ::QuantumObject{ObjType};
    base::Int = 0,
    tol::Real = 1e-15,
) where {ObjType<:Union{KetQuantumObject,OperatorQuantumObject}}
    T = eltype(ρ)
    vals = eigenenergies(ket2dm(ρ))
    indexes = findall(x -> abs(x) > tol, vals)
    length(indexes) == 0 && return zero(real(T))
    nzvals = vals[indexes]
    logvals = base != 0 ? log.(base, Complex.(nzvals)) : log.(Complex.(nzvals))
    return -real(mapreduce(*, +, nzvals, logvals))
end

"""
    entanglement(QO::QuantumObject, sel::Union{Int,AbstractVector{Int},Tuple})

Calculates the entanglement by doing the partial trace of `QO`, selecting only the dimensions
with the indices contained in the `sel` vector, and then using the Von Neumann entropy [`entropy_vn`](@ref).
"""
function entanglement(
    QO::QuantumObject{OpType},
    sel::Union{AbstractVector{Int},Tuple},
) where {OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    ψ = normalize(QO)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end
entanglement(QO::QuantumObject, sel::Int) = entanglement(QO, (sel,))
