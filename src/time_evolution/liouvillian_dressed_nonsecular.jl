export liouvillian_dressed_nonsecular

@doc raw"""
    liouvillian_dressed_nonsecular(
        H::QuantumObject{Operator},
        fields::Vector,
        T_list::Vector{<:Real};
        N_trunc::Union{Int,Nothing}=nothing,
        tol::Real=1e-12,
        σ_filter::Union{Nothing,Real}=nothing,
    )

Build the generalized Liouvillian for a system coupled to multiple bosonic baths in the ultrastrong-coupling regime. The Hamiltonian `H` is diagonalized, the system-bath operators in `fields` are projected into the eigenbasis, and thermal jump operators are assembled for each temperature in `T_list`.

# Arguments
- `H::QuantumObject{Operator}`: System Hamiltonian.
- `fields::Vector`: Coupling operators that mediate the interaction with each bath; must match the length of `T_list`.
- `T_list::Vector{<:Real}`: Bath temperatures ordered consistently with `fields`.

# Keyword Arguments
- `N_trunc::Union{Int,Nothing}`: If provided, truncate the eigenbasis to the first `N_trunc` levels; otherwise use the full dimension of `H`.
- `tol::Real`: Tolerance passed to sparsification utilities when constructing the filters and dissipators.
- `σ_filter::Union{Nothing,Real}`: Width of the Gaussian frequency filter. If `nothing`, a heuristic value proportional to the coupling strengths is used.

# Returns
- `E`: Eigenenergies of `H` (truncated if `N_trunc` is provided).
- `U`: Eigenvectors of `H` as a [`QuantumObject`](@ref) mapping to the truncated basis.
- `L`: Generalized Liouvillian [`SuperOperator`](@ref) including the frequency-filtered dissipators.

# References
- [Settineri2018](@cite)
"""
function liouvillian_dressed_nonsecular(
        H::QuantumObject{Operator},
        fields::Vector,
        T_list::Vector{<:Real};
        N_trunc::Union{Int, Nothing} = nothing,
        tol::Real = 1.0e-12,
        σ_filter::Union{Nothing, Real} = nothing,
    )
    (length(fields) == length(T_list)) || throw(DimensionMismatch("The number of fields and T_list must be the same."))

    dims = isnothing(N_trunc) ? H.dimensions : ProductDimensions(N_trunc)
    final_size = get_hilbert_size(dims)[1]
    # U is a non-square transformation matrix from original basis to truncated eigenbasis
    final_dims = isnothing(N_trunc) ? H.dimensions : ProductDimensions(H.dimensions.to, (HilbertSpace(N_trunc),))
    result = eigen(H)
    E = real.(result.values[1:final_size])
    U = QuantumObject(result.vectors[:, 1:final_size], result.type, final_dims)

    H_d = QuantumObject(Diagonal(complex(E)), type = Operator(), dims = dims)

    Ω = E' .- E
    Ωp = triu(to_sparse(Ω, tol), 1)

    # Filter width
    σ = isnothing(σ_filter) ? 500 * maximum([norm(field) / length(field) for field in fields]) : σ_filter

    L = liouvillian(H_d) + sum(eachindex(fields)) do i
        # The operator that couples the system to the bath in the eigenbasis
        X_op = to_sparse((U' * fields[i] * U).data, tol)
        if ishermitian(fields[i])
            X_op = (X_op + X_op') / 2 # Make sure it's hermitian
        end

        # Ohmic reservoir
        N_th = n_thermal.(Ωp, T_list[i])
        Sp₀ = QuantumObject(triu(X_op, 1), type = Operator(), dims = dims)
        Sp₁ = QuantumObject(droptol!((@. Ωp * N_th * Sp₀.data), tol), type = Operator(), dims = dims)
        Sp₂ = QuantumObject(droptol!((@. Ωp * (1 + N_th) * Sp₀.data), tol), type = Operator(), dims = dims)
        # S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

        # Build the dissipator contribution with filtered spre, spost, sprepost to avoid large intermediate kron allocations.
        D₁ = (
            _filtered_sprepost(Sp₁', Sp₀, E, σ, tol) +
                _filtered_sprepost(Sp₀', Sp₁, E, σ, tol) -
                _filtered_spre(Sp₀ * Sp₁', E, σ, tol) -
                _filtered_spost(Sp₁ * Sp₀', E, σ, tol)
        ) / 2
        D₂ = (
            _filtered_sprepost(Sp₂, Sp₀', E, σ, tol) +
                _filtered_sprepost(Sp₀, Sp₂', E, σ, tol) -
                _filtered_spre(Sp₀' * Sp₂, E, σ, tol) -
                _filtered_spost(Sp₂' * Sp₀, E, σ, tol)
        ) / 2

        D₁ + D₂
    end

    settings.auto_tidyup && tidyup!(L)

    return E, U, L
end

function _filtered_sprepost(
        A::AbstractQuantumObject{Operator},
        B::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
    )
    isinf(σ) && return sprepost(A, B)

    check_dimensions(A, B)
    data = _filtered_kron(transpose(B.data), A.data, E, σ, tol)
    return promote_op_type(A, B)(data, SuperOperator(), A.dimensions)
end

function _filtered_spre(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
    )
    isinf(σ) && return spre(A)

    T = eltype(A)
    data = _filtered_kron(Eye{T}(size(A, 1)), A.data, E, σ, tol)
    return get_typename_wrapper(A)(data, SuperOperator(), A.dimensions)
end

function _filtered_spost(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
    )
    isinf(σ) && return spost(A)

    T = eltype(A)
    data = _filtered_kron(transpose(A.data), Eye{T}(size(A, 1)), E, σ, tol)
    return get_typename_wrapper(A)(data, SuperOperator(), A.dimensions)
end

function _filtered_kron(A, B, E, σ, tol)
    N = length(E)
    (size(A, 1) == N && size(A, 2) == N && size(B, 1) == N && size(B, 2) == N) ||
        throw(DimensionMismatch("Matrix sizes do not match energy list; expected $(N)×$(N) matrices."))

    I_A, J_A, V_A = _findnz(A)
    I_B, J_B, V_B = _findnz(B)

    nnzA = length(V_A)
    nnzB = length(V_B)

    T = Base.promote_eltype(V_A, V_B, E, σ)
    I_out = Int[]
    J_out = Int[]
    V_out = Vector{T}()

    (nnzA == 0 || nnzB == 0) && return sparse(I_out, J_out, V_out, N^2, N^2)

    sizehint!(I_out, max(nnzA, nnzB))
    sizehint!(J_out, max(nnzA, nnzB))
    sizehint!(V_out, max(nnzA, nnzB))

    @inbounds for idxA in eachindex(V_A)
        iA = I_A[idxA]
        jA = J_A[idxA]
        vA = V_A[idxA]
        for idxB in eachindex(V_B)
            iT = I_B[idxB]
            jT = J_B[idxB]
            vB = V_B[idxB]
            Ωdiff = (E[iA] - E[jA]) - (E[iT] - E[jT])
            v = vA * vB * gaussian(Ωdiff, 0, σ)
            if abs(v) > tol
                push!(I_out, iA + (iT - 1) * N)
                push!(J_out, jA + (jT - 1) * N)
                push!(V_out, v)
            end
        end
    end

    return sparse(I_out, J_out, V_out, N^2, N^2)
end
function _filtered_kron(A::Transpose{T, <:Adjoint}, B, E, σ, tol) where {T}
    return _filtered_kron(conj(parent(parent(A))), B, E, σ, tol)
end

_findnz(A::AbstractSparseMatrix) = findnz(A)
function _findnz(A::Diagonal{T}) where {T}
    V = diag(A)
    I = similar(V, Int)
    I .= 1:length(V)
    J = I
    return I, J, V
end
function _findnz(A::Transpose{T, <:AbstractMatrix}) where {T}
    I, J, V = findnz(parent(A))

    return J, I, V
end
function _findnz(A::Adjoint{T, <:AbstractMatrix}) where {T}
    I, J, V = findnz(parent(A))

    return J, I, conj.(V)
end
