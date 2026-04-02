#=
Entropy related functions and some entanglement measures.
=#

export entropy_vn, entropy_relative, entropy_linear, entropy_mutual, entropy_conditional
export entanglement, concurrence

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

Quantum Object:   type=Ket()   dims=([2], [1])   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im

julia> entropy_vn(ψ, base=2)
-0.0
```

Mixed state:
```jldoctest
julia> ρ = maximally_mixed_dm(2)

Quantum Object:   type=Operator()   dims=([2], [2])   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im

julia> entropy_vn(ρ, base=2)
1.0
```
"""
function entropy_vn(ρ::QuantumObject{ObjType}; base::Int = 0, tol::Real = 1.0e-15) where {ObjType <: Union{Ket, Operator}}
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
        tol::Real = 1.0e-15,
    ) where {ObjType1 <: Union{Ket, Operator}, ObjType2 <: Union{Ket, Operator}}
    ρ_dm = ket2dm(ρ)
    σ_dm = ket2dm(σ)
    check_dimensions(ρ_dm, σ_dm)

    # the logic of this code follows the detail given in the reference of the docstring
    # consider the eigen decompositions:
    #   ρ = Σ_i p_i |i⟩⟨i|
    #   σ = Σ_j q_j |j⟩⟨j|
    ρ_result = eigenstates(ρ_dm)
    σ_result = eigenstates(σ_dm)

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
    return max(0.0, dot(p_vals, log_p) - dot(p, P, log_q)) # use 0.0 to make sure it always return value in Float-type
end

@doc raw"""
    entropy_linear(ρ::QuantumObject)

Calculates the quantum linear entropy ``S_L = 1 - \textrm{Tr} \left[ \hat{\rho}^2 \right]``, where ``\hat{\rho}`` is the density matrix of the system.

Note that `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref).
"""
entropy_linear(ρ::QuantumObject{ObjType}) where {ObjType <: Union{Ket, Operator}} = 1.0 - purity(ρ) # use 1.0 to make sure it always return value in Float-type

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
        ρAB::QuantumObject{ObjType, <:Dimensions{<:TensorSpace{N}, <:TensorSpace{N}}}, # the dimensions to == from, and should both be TensorSpace
        selA::Union{Int, AbstractVector{Int}, Tuple},
        selB::Union{Int, AbstractVector{Int}, Tuple};
        kwargs...,
    ) where {ObjType <: Union{Ket, Operator}, N}
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
    ρAB::QuantumObject{ObjType, <:Dimensions{N, N}},
    selB::Union{Int, AbstractVector{Int}, Tuple};
    kwargs...,
) where {ObjType <: Union{Ket, Operator}, N} = entropy_vn(ρAB; kwargs...) - entropy_vn(ptrace(ρAB, selB); kwargs...)

@doc raw"""
    entanglement(ρ::QuantumObject, sel; kwargs...)

Calculates the [entanglement entropy](https://en.wikipedia.org/wiki/Entropy_of_entanglement) by doing the partial trace of `ρ`, selecting only the dimensions with the indices contained in the `sel` vector, and then use the Von Neumann entropy [`entropy_vn`](@ref).

# Notes

- `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref). But should be a pure state.
- `sel` specifies the indices of the remaining sub-system. See also [`ptrace`](@ref).
- `kwargs` are the keyword arguments for calculating Von Neumann entropy. See also [`entropy_vn`](@ref).
"""
function entanglement(
        ρ::QuantumObject{OpType},
        sel::Union{Int, AbstractVector{Int}, Tuple},
        kwargs...,
    ) where {OpType <: Union{Ket, Operator}}
    p = purity(ρ)
    isapprox(p, 1; atol = 1.0e-2) || throw(
        ArgumentError(
            "The entanglement entropy only works for normalized pure state, the purity of the given state: $(p) ≉ 1",
        ),
    )

    ρ_tr = ptrace(ρ, sel)
    val = entropy_vn(ρ_tr; kwargs...)
    return max(0.0, val)  # use 0.0 to make sure it always return value in Float-type
end

@doc raw"""
    concurrence(ρ::QuantumObject)

Calculate the [concurrence](https://en.wikipedia.org/wiki/Concurrence_(quantum_computing)) for a two-qubit state.

# Notes

- `ρ` can be either a [`Ket`](@ref) or an [`Operator`](@ref).

# References

- [Hill-Wootters1997](@citet)
"""
function concurrence(ρ::QuantumObject{OpType}) where {OpType <: Union{Ket, Operator}}
    two_qubit_dims = TensorSpace(Space(2), Space(2))
    is_two_qubit = (isket(ρ) || isendomorphic(ρ.dimensions)) && ρ.dimensions.to == two_qubit_dims
    is_two_qubit || throw(
        ArgumentError(
            "The `concurrence` only works for a two-qubit state, invalid dims = $(_get_dims_string(ρ.dimensions)).",
        ),
    )

    ρ_mat = ket2dm(ρ).data
    σyσy = tensor(sigmay(), sigmay()).data
    ρ_tilde = σyσy * conj(ρ_mat) * σyσy

    # We use the alternative way to calculate concurrence (more efficient):
    # calculate the square root of each eigenvalue (in decreasing order) of the non-Hermitian matrix ρ * ρ̃.
    # Note: we use abs() to avoid problems with sqrt for very small negative numbers due to numerical precision.
    λ = sqrt.(abs.(real(eigvals(ρ_mat * ρ_tilde; sortby = x -> -real(x)))))

    return max(zero(λ[1]), λ[1] - λ[2] - λ[3] - λ[4])
end
