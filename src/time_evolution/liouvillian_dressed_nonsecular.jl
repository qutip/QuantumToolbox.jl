export liouvillian_dressed_nonsecular

@doc raw"""
    liouvillian_dressed_nonsecular(
        H::QuantumObject{Operator},
        fields::Vector,
        T_list::Vector{<:Real};
        N_trunc::Union{Int,Nothing}=nothing,
        tol::Real=1e-12,
        σ_filter::Union{Nothing,Real}=nothing,
        matrix_form::Union{Bool,Val}=Val(false),
    )

Build the generalized Liouvillian for a system coupled to multiple bosonic baths when the system subparts are highly coupled. The Hamiltonian `H` is diagonalized, the system-bath operators in `fields` are projected into the eigenbasis, and thermal jump operators are assembled for each temperature in `T_list`.

# Arguments
- `H::QuantumObject{Operator}`: System Hamiltonian.
- `fields::Vector`: Coupling operators that mediate the interaction with each bath; must match the length of `T_list`.
- `T_list::Vector{<:Real}`: Bath temperatures ordered consistently with `fields`.

# Keyword Arguments
- `N_trunc::Union{Int,Nothing}`: If provided, truncate the eigenbasis to the first `N_trunc` levels; otherwise use the full dimension of `H`.
- `tol::Real`: Tolerance passed to sparsification utilities when constructing the filters and dissipators.
- `σ_filter::Union{Nothing,Real}`: Width of the Gaussian frequency filter. If `nothing`, filtering is disabled.
- `matrix_form::Union{Bool,Val}`: If `Val(true)`, return the Liouvillian in matrix form (`SuperOperatorMatrixForm`); otherwise return the vectorized `SuperOperator`.

# Returns
- `E`: Eigenenergies of `H` (truncated if `N_trunc` is provided).
- `U`: Eigenvectors of `H` as a [`QuantumObject`](@ref) mapping to the truncated basis.
- `L`: Generalized Liouvillian [`SuperOperator`](@ref) including the frequency-filtered dissipators.

# References
- [Settineri2018](@cite)
"""
function liouvillian_dressed_nonsecular(
        H::QuantumObject{Operator},
        fields::Base.AbstractVecOrTuple,
        T_list::Base.AbstractVecOrTuple{<:Real};
        N_trunc::Union{Int, Nothing} = nothing,
        tol::Real = 1.0e-12,
        σ_filter::Union{Nothing, Real} = nothing,
        matrix_form::Union{Bool, Val} = Val(false),
    )
    (length(fields) == length(T_list)) || throw(DimensionMismatch("The number of fields and T_list must be the same."))

    dims = isnothing(N_trunc) ? H.dimensions : Dimensions(N_trunc)
    final_size = get_size(dims)[1]
    # U is a non-square transformation matrix from original basis to truncated eigenbasis
    final_dims = isnothing(N_trunc) ? H.dimensions : Dimensions(H.dimensions.to, Space(N_trunc))
    result = eigen(H)
    E = real.(result.values[1:final_size])
    U = QuantumObject(result.vectors[:, 1:final_size], result.type, final_dims)

    H_d = QuantumObject(Diagonal(complex(E)), type = Operator(), dims = dims)

    Ω = E' .- E
    Ωp = triu(to_sparse(Ω, tol), 1)

    # mapreduce doesn't work with tuples, so we do map and reduce separately
    S_spre, S_spost, L_sprepost = reduce(
        .+, map(T_list, fields) do T, field
            # The operator that couples the system to the bath in the eigenbasis
            X_op = to_sparse((U' * field * U).data, tol)
            if ishermitian(field)
                X_op = (X_op + X_op') / 2 # Make sure it's hermitian
            end

            # Ohmic reservoir
            N_th = n_thermal.(Ωp, T)
            Sp₀ = QuantumObject(triu(X_op, 1), type = Operator(), dims = dims)
            Sp₁ = QuantumObject(droptol!((@. Ωp * N_th * Sp₀.data), tol), type = Operator(), dims = dims)
            Sp₂ = QuantumObject(droptol!((@. Ωp * (1 + N_th) * Sp₀.data), tol), type = Operator(), dims = dims)
            # S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

            S_spre_i = - Sp₀ * Sp₁' / 2 - Sp₀' * Sp₂ / 2
            S_spost_i = - Sp₁ * Sp₀' / 2 - Sp₂' * Sp₀ / 2
            S_sprepost_i = _filtered_sprepost(Sp₁', Sp₀ / 2, E, σ_filter, tol, matrix_form) +
                _filtered_sprepost(Sp₀', Sp₁ / 2, E, σ_filter, tol, matrix_form) +
                _filtered_sprepost(Sp₂, Sp₀' / 2, E, σ_filter, tol, matrix_form) +
                _filtered_sprepost(Sp₀, Sp₂' / 2, E, σ_filter, tol, matrix_form)

            return S_spre_i, S_spost_i, S_sprepost_i
        end
    )

    L_spre = _filtered_spre(-im * H_d + S_spre, E, σ_filter, tol, matrix_form)
    L_spost = _filtered_spost(im * H_d + S_spost, E, σ_filter, tol, matrix_form)

    L = L_spre + L_spost + L_sprepost

    settings.auto_tidyup && (L isa QuantumObject) && tidyup!(L) # tidyup! only supports QuantumObject

    return E, U, L
end

function _filtered_sprepost(
        A::AbstractQuantumObject{Operator},
        B::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Nothing,
        tol::Real,
        matrix_form,
    )
    return sprepost(A, B; matrix_form = matrix_form)
end

function _filtered_sprepost(
        A::AbstractQuantumObject{Operator},
        B::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
        matrix_form,
    )
    getVal(matrix_form) &&
        throw(ArgumentError("Filtered matrix-form Liouvillian is not implemented yet. Use `σ_filter = nothing` for unfiltered matrix-form output."))

    check_dimensions(A, B)
    data = _filtered_kron(transpose(B.data), A.data, E, σ, tol)
    return promote_op_type(A, B)(data, SuperOperator(), LiouvilleSpace(A.dimensions))
end

function _filtered_spre(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Nothing,
        tol::Real,
        matrix_form,
    )
    return spre(A; matrix_form = matrix_form)
end

function _filtered_spre(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
        matrix_form,
    )
    getVal(matrix_form) &&
        throw(ArgumentError("Filtered matrix-form Liouvillian is not implemented yet. Use `σ_filter = nothing` for unfiltered matrix-form output."))

    T = eltype(A)
    data = _filtered_kron(Eye{T}(size(A, 1)), A.data, E, σ, tol)
    return get_typename_wrapper(A)(data, SuperOperator(), LiouvilleSpace(A.dimensions))
end

function _filtered_spost(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Nothing,
        tol::Real,
        matrix_form,
    )
    return spost(A; matrix_form = matrix_form)
end

function _filtered_spost(
        A::AbstractQuantumObject{Operator},
        E::AbstractVector,
        σ::Real,
        tol::Real,
        matrix_form,
    )
    getVal(matrix_form) &&
        throw(ArgumentError("Filtered matrix-form Liouvillian is not implemented yet. Use `σ_filter = nothing` for unfiltered matrix-form output."))

    T = eltype(A)
    data = _filtered_kron(transpose(A.data), Eye{T}(size(A, 1)), E, σ, tol)
    return get_typename_wrapper(A)(data, SuperOperator(), LiouvilleSpace(A.dimensions))
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

    sizehint!(I_out, max(nnzA, nnzB, div(nnzA * nnzB, 10)))
    sizehint!(J_out, max(nnzA, nnzB, div(nnzA * nnzB, 10)))
    sizehint!(V_out, max(nnzA, nnzB, div(nnzA * nnzB, 10)))

    @inbounds for idxA in eachindex(V_A)
        iA = I_A[idxA]
        jA = J_A[idxA]
        vA = V_A[idxA]

        for idxB in eachindex(V_B)
            iB = I_B[idxB]
            jB = J_B[idxB]
            vB = V_B[idxB]

            Ωdiff = (E[iA] - E[jA]) - (E[iB] - E[jB])
            v = vA * vB * gaussian(Ωdiff, 0, σ)

            if abs(v) > tol
                # Corrected mapping for A ⊗ B
                push!(I_out, iB + (iA - 1) * N)
                push!(J_out, jB + (jA - 1) * N)
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
