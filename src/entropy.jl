#=
Entropy related functions.
=#

export entropy_vn, entropy_linear, entropy_mutual, entropy_conditional
export entanglement

@doc raw"""
    entropy_vn(ρ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy) ``S = - \textrm{Tr} \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right]``, where ``\hat{\rho}`` is the density matrix of the system.

# Notes

- `ρ` is the quantum state, can be either a [`Ket`](@ref) or an [`Operator`](@ref).
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

@doc raw"""
    entropy_linear(ρ::QuantumObject)

Calculates the linear entropy ``S_L = 1 - \textrm{Tr} \left[ \hat{\rho}^2 \right]``, where ``\hat{\rho}`` is the density matrix of the system.

Note that `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
"""
entropy_linear(ρ::QuantumObject{ObjType}) where {ObjType<:Union{KetQuantumObject,OperatorQuantumObject}} =
    1.0 - purity(ρ) # use 1.0 to make sure it always return value in Float-type

@doc raw"""
    entropy_mutual(ρAB::QuantumObject, selA, selB; kwargs...)

Calculates the mutual information ``I(A:B) = S(\hat{\rho}_A) + S(\hat{\rho}_B) - S(\hat{\rho}_{AB})`` between subsystems ``A`` and ``B``.

Here, ``S`` is the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy), ``\hat{\rho}_{AB}`` is the density matrix of the entire system, ``\hat{\rho}_A = \textrm{Tr}_B \left[ \hat{\rho}_{AB} \right]``, ``\hat{\rho}_B = \textrm{Tr}_A \left[ \hat{\rho}_{AB} \right]``.

# Notes

- `ρAB` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `selA` specifies the indices of the sub-system `A` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `selB` specifies the indices of the sub-system `B` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
function entropy_mutual(
    ρAB::QuantumObject{ObjType,<:AbstractDimensions{N}},
    selA::AType,
    selB::BType;
    kwargs...,
) where {
    ObjType<:Union{KetQuantumObject,OperatorQuantumObject},
    N,
    AType<:Union{Int,AbstractVector{Int},Tuple},
    BType<:Union{Int,AbstractVector{Int},Tuple},
}
    # check if selA and selB matches the dimensions of ρAB
    sel_A_B = vcat(selA, selB)
    (length(sel_A_B) != N) && ArgumentError(
        "The indices in `selA = $(selA)` and `selB = $(selB)` does not match the given QuantumObject which has $N sub-systems",
    )
    allunique(sel_A_B) || throw(ArgumentError("Duplicate selection indices in `selA = $(selA)` and `selB = $(selB)`"))

    ρA = ptrace(ρAB, selA)
    ρB = ptrace(ρAB, selB)
    return entropy_vn(ρA; kwargs...) + entropy_vn(ρB; kwargs...) - entropy_vn(ρAB; kwargs...)
end

@doc raw"""
    entropy_conditional(ρAB::QuantumObject, selB; kwargs...)

Calculates the conditional entropy with respect to sub-system ``B``: ``S(A|B) = S(\hat{\rho}_{AB}) - S(\hat{\rho}_{B})``.

Here, ``S`` is the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy), ``\hat{\rho}_{AB}`` is the density matrix of the entire system, and ``\hat{\rho}_B = \textrm{Tr}_A \left[ \hat{\rho}_{AB} \right]``.

# Notes

- `ρAB` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `selB` specifies the indices of the sub-system `B` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
entropy_conditional(
    ρAB::QuantumObject{ObjType,<:AbstractDimensions{N}},
    selB::BType;
    kwargs...,
) where {
    ObjType<:Union{KetQuantumObject,OperatorQuantumObject},
    N,
    BType<:Union{Int,AbstractVector{Int},Tuple},
} = entropy_vn(ρAB; kwargs...) - entropy_vn(ptrace(ρAB, selB); kwargs...)

"""
    entanglement(QO::QuantumObject, sel::Union{Int,AbstractVector{Int},Tuple})

Calculates the entanglement by doing the partial trace of `QO`, selecting only the dimensions with the indices contained in the `sel` vector, and then using the Von Neumann entropy [`entropy_vn`](@ref).
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
