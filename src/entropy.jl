#=
Entropy related functions.
=#

export entropy_vn, entropy_relative, entropy_linear, entropy_mutual, entropy_conditional
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
    entropy_relative(ρ::QuantumObject, σ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [quantum relative entropy](https://en.wikipedia.org/wiki/Quantum_relative_entropy) of ``\hat{\rho}`` with respect to ``\hat{\sigma}``: ``D(\hat{\rho}||\hat{\sigma}) = \textrm{Tr} \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right] - \textrm{Tr} \left[ \hat{\rho} \log \left( \hat{\sigma} \right) \right]``.

# Notes

- `ρ` is a quantum state, can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `σ` is a quantum state, can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `base` specifies the base of the logarithm to use, and when using the default value `0`, the natural logarithm is used.
- `tol` describes the absolute tolerance for detecting the zero-valued eigenvalues of the density matrix ``\hat{\rho}``.

# References
- [Nielsen-Chuang2011; section 11.3.1, page 511](@citet)
"""
function entropy_relative(
    ρ::QuantumObject{ObjType1},
    σ::QuantumObject{ObjType2};
    base::Int = 0,
    tol::Real = 1e-15,
) where {
    ObjType1<:Union{KetQuantumObject,OperatorQuantumObject},
    ObjType2<:Union{KetQuantumObject,OperatorQuantumObject},
}
    check_dimensions(ρ, σ)

    # the logic of this code follows the detail given in the reference of the docstring
    # consider the eigen decompositions:
    #   ρ = Σ_i p_i |i⟩⟨i|
    #   σ = Σ_j q_j |j⟩⟨j|
    ρ_result = eigenstates(ket2dm(ρ))
    σ_result = eigenstates(ket2dm(σ))

    # make sure all p_i and q_j are real
    any(p_i -> imag(p_i) >= tol, ρ_result.values) && error("Input `ρ` has non-real eigenvalues.")
    any(q_j -> imag(q_j) >= tol, σ_result.values) && error("Input `σ` has non-real eigenvalues.")
    p = real(ρ_result.values)
    q = real(σ_result.values)
    Uρ = ρ_result.vectors
    Uσ = σ_result.vectors

    # create P_ij matrix (all elements should be real)
    P = abs2.(Uρ' * Uσ) # this equals to ⟨i|j⟩⟨j|i⟩

    # return +∞ if kernel of σ overlaps with support of ρ, i.e., supp(p) ⊆ supp(q)
    # That is, if σ is not full rank, S(ρ||σ) = +∞
    # note that, one special case is that S(ρ||σ) = 0 (if ρ == σ)
    ((transpose(p .>= tol) * (P .>= tol) * (q .< tol)) == 0) || return Inf

    # Avoid -∞ from log(0), these terms will be multiplied by zero later anyway
    replace!(q_j -> abs(q_j) < tol ? 1 : q_j, q)
    p_vals = filter(p_i -> abs(p_i) >= tol, p)

    if base == 0
        log_p = log.(p_vals)
        log_q = log.(q)
    else
        log_p = log.(base, p_vals)
        log_q = log.(base, q)
    end

    # the relative entropy is guaranteed to be ≥ 0
    # so we calculate the value to 0 to avoid small violations of the lower bound.
    return max(0.0, dot(p_vals, log_p) - dot(p, P, log_q))
end

@doc raw"""
    entropy_linear(ρ::QuantumObject)

Calculates the quantum linear entropy ``S_L = 1 - \textrm{Tr} \left[ \hat{\rho}^2 \right]``, where ``\hat{\rho}`` is the density matrix of the system.

Note that `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
"""
entropy_linear(ρ::QuantumObject{ObjType}) where {ObjType<:Union{KetQuantumObject,OperatorQuantumObject}} =
    1.0 - purity(ρ) # use 1.0 to make sure it always return value in Float-type

@doc raw"""
    entropy_mutual(ρAB::QuantumObject, selA, selB; kwargs...)

Calculates the [quantum mutual information](https://en.wikipedia.org/wiki/Quantum_mutual_information) ``I(A:B) = S(\hat{\rho}_A) + S(\hat{\rho}_B) - S(\hat{\rho}_{AB})`` between subsystems ``A`` and ``B``.

Here, ``S`` is the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy), ``\hat{\rho}_{AB}`` is the density matrix of the entire system, ``\hat{\rho}_A = \textrm{Tr}_B \left[ \hat{\rho}_{AB} \right]``, ``\hat{\rho}_B = \textrm{Tr}_A \left[ \hat{\rho}_{AB} \right]``.

# Notes

- `ρAB` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `selA` specifies the indices of the sub-system `A` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `selB` specifies the indices of the sub-system `B` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
function entropy_mutual(
    ρAB::QuantumObject{ObjType,<:AbstractDimensions{N}},
    selA::Union{Int,AbstractVector{Int},Tuple},
    selB::Union{Int,AbstractVector{Int},Tuple};
    kwargs...,
) where {ObjType<:Union{KetQuantumObject,OperatorQuantumObject},N}
    # check if selA and selB matches the dimensions of ρAB
    sel_A_B = (selA..., selB...)
    (length(sel_A_B) != N) && throw(
        ArgumentError(
            "The indices in `selA = $(selA)` and `selB = $(selB)` does not match the given QuantumObject which has $N sub-systems",
        ),
    )
    allunique(sel_A_B) || throw(ArgumentError("Duplicate selection indices in `selA = $(selA)` and `selB = $(selB)`"))

    ρA = ptrace(ρAB, selA)
    ρB = ptrace(ρAB, selB)
    return entropy_vn(ρA; kwargs...) + entropy_vn(ρB; kwargs...) - entropy_vn(ρAB; kwargs...)
end

@doc raw"""
    entropy_conditional(ρAB::QuantumObject, selB; kwargs...)

Calculates the [conditional quantum entropy](https://en.wikipedia.org/wiki/Conditional_quantum_entropy) with respect to sub-system ``B``: ``S(A|B) = S(\hat{\rho}_{AB}) - S(\hat{\rho}_{B})``.

Here, ``S`` is the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy), ``\hat{\rho}_{AB}`` is the density matrix of the entire system, and ``\hat{\rho}_B = \textrm{Tr}_A \left[ \hat{\rho}_{AB} \right]``.

# Notes

- `ρAB` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `selB` specifies the indices of the sub-system `B` in `ρAB.dimensions`. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
entropy_conditional(
    ρAB::QuantumObject{ObjType,<:AbstractDimensions{N}},
    selB::Union{Int,AbstractVector{Int},Tuple};
    kwargs...,
) where {ObjType<:Union{KetQuantumObject,OperatorQuantumObject},N} =
    entropy_vn(ρAB; kwargs...) - entropy_vn(ptrace(ρAB, selB); kwargs...)

@doc raw"""
    entanglement(ρ::QuantumObject, sel; kwargs...)

Calculates the [entanglement entropy](https://en.wikipedia.org/wiki/Entropy_of_entanglement) by doing the partial trace of `ρ`, selecting only the dimensions with the indices contained in the `sel` vector, and then use the Von Neumann entropy [`entropy_vn`](@ref).

# Notes

- `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
- `sel` specifies the indices of the remaining sub-system. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
function entanglement(
    ρ::QuantumObject{OpType},
    sel::Union{Int,AbstractVector{Int},Tuple},
    kwargs...,
) where {OpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    _ρ = normalize(ρ)
    ρ_tr = ptrace(_ρ, sel)
    val = entropy_vn(ρ_tr; kwargs...)
    return (val > 0) * val
end
