export OperatorSum, TimeDependentOperatorSum
export TimeEvolutionSol, TimeEvolutionMCSol

export liouvillian, liouvillian_floquet, liouvillian_generalized

struct TimeEvolutionSol{TT<:Vector{<:Real},TS<:AbstractVector,TE<:Matrix{ComplexF64}}
    times::TT
    states::TS
    expect::TE
end

struct TimeEvolutionMCSol{
    TT<:Vector{<:Vector{<:Real}},
    TS<:AbstractVector,
    TE<:Matrix{ComplexF64},
    TEA<:Array{ComplexF64,3},
    TJT<:Vector{<:Vector{<:Real}},
    TJW<:Vector{<:Vector{<:Integer}},
}
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    jump_times::TJT
    jump_which::TJW
end

abstract type LindbladJumpCallbackType end

struct ContinuousLindbladJumpCallback <: LindbladJumpCallbackType
    interp_points::Int
end

struct DiscreteLindbladJumpCallback <: LindbladJumpCallbackType end

ContinuousLindbladJumpCallback(; interp_points::Int = 10) = ContinuousLindbladJumpCallback(interp_points)

## Sum of operators

struct OperatorSum{CT<:Vector{<:Number},OT<:Vector{<:QuantumObject}} <: AbstractQuantumObject
    coefficients::CT
    operators::OT
    function OperatorSum(coefficients::CT, operators::OT) where {CT<:Vector{<:Number},OT<:Vector{<:QuantumObject}}
        length(coefficients) == length(operators) ||
            throw(DimensionMismatch("The number of coefficients must be the same as the number of operators."))
        # Check if all the operators have the same dimensions
        dims = operators[1].dims
        optype = operators[1].type
        mapreduce(x -> x.dims == dims && x.type == optype, &, operators) ||
            throw(DimensionMismatch("All the operators must have the same dimensions."))
        T = promote_type(
            mapreduce(x -> eltype(x.data), promote_type, operators),
            mapreduce(eltype, promote_type, coefficients),
        )
        coefficients2 = T.(coefficients)
        return new{Vector{T},OT}(coefficients2, operators)
    end
end

Base.size(A::OperatorSum) = size(A.operators[1])
Base.size(A::OperatorSum, inds...) = size(A.operators[1], inds...)
Base.length(A::OperatorSum) = length(A.operators[1])
Base.copy(A::OperatorSum) = OperatorSum(copy(A.coefficients), copy(A.operators))
Base.deepcopy(A::OperatorSum) = OperatorSum(deepcopy(A.coefficients), deepcopy(A.operators))

function update_coefficients!(A::OperatorSum, coefficients)
    length(A.coefficients) == length(coefficients) ||
        throw(DimensionMismatch("The number of coefficients must be the same as the number of operators."))
    return A.coefficients .= coefficients
end

@inline function LinearAlgebra.mul!(y::AbstractVector{T}, A::OperatorSum, x::AbstractVector, α, β) where {T}
    # Note that β is applied only to the first term
    mul!(y, A.operators[1], x, α * A.coefficients[1], β)
    @inbounds for i in 2:length(A.operators)
        A.coefficients[i] == 0 && continue
        mul!(y, A.operators[i], x, α * A.coefficients[i], 1)
    end
    return y
end

struct TimeDependentOperatorSum{CFT,OST<:OperatorSum}
    coefficient_functions::CFT
    operator_sum::OST
end

function TimeDependentOperatorSum(
    coefficient_functions,
    operators::Vector{<:QuantumObject};
    params = nothing,
    init_time = 0.0,
)
    # promote the type of the coefficients and the operators. Remember that the coefficient_functions si a vector of functions and the operators is a vector of QuantumObjects
    coefficients = [f(init_time, params) for f in coefficient_functions]
    operator_sum = OperatorSum(coefficients, operators)
    return TimeDependentOperatorSum(coefficient_functions, operator_sum)
end

Base.size(A::TimeDependentOperatorSum) = size(A.operator_sum)
Base.size(A::TimeDependentOperatorSum, inds...) = size(A.operator_sum, inds...)
Base.length(A::TimeDependentOperatorSum) = length(A.operator_sum)

function update_coefficients!(A::TimeDependentOperatorSum, t, params)
    @inbounds @simd for i in 1:length(A.coefficient_functions)
        A.operator_sum.coefficients[i] = A.coefficient_functions[i](t, params)
    end
end

(A::TimeDependentOperatorSum)(t, params) = (update_coefficients!(A, t, params); A)

@inline function LinearAlgebra.mul!(y::AbstractVector, A::TimeDependentOperatorSum, x::AbstractVector, α, β)
    return mul!(y, A.operator_sum, x, α, β)
end

#######################################

### LIOUVILLIAN ###
@doc raw"""
    liouvillian(H::QuantumObject, c_ops::AbstractVector, Id_cache=I(prod(H.dims))

Construct the Liouvillian superoperator for a system Hamiltonian and a set of collapse operators:
``\mathcal{L} \cdot = -i[\hat{H}, \cdot] + \sum_i \mathcal{D}[\hat{O}_i] \cdot``,
where ``\mathcal{D}[\hat{O}_i] \cdot = \hat{O}_i \cdot \hat{O}_i^\dagger - \frac{1}{2} \hat{O}_i^\dagger \hat{O}_i \cdot - \frac{1}{2} \cdot \hat{O}_i^\dagger \hat{O}_i``.

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when
the same function is applied multiple times with a known Hilbert space dimension.
"""
function liouvillian(
    H::QuantumObject{MT1,OpType1},
    c_ops::Vector{QuantumObject{MT2,OpType2}} = Vector{QuantumObject{MT1,OpType1}}([]),
    Id_cache = I(prod(H.dims)),
) where {
    MT1<:AbstractMatrix,
    MT2<:AbstractMatrix,
    OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
}
    L = liouvillian(H, Id_cache)
    for c_op in c_ops
        L += lindblad_dissipator(c_op, Id_cache)
    end
    return L
end

liouvillian(
    H::QuantumObject{MT1,OperatorQuantumObject},
    Id_cache::Diagonal = I(prod(H.dims)),
) where {MT1<:AbstractMatrix} = -1im * (spre(H, Id_cache) - spost(H, Id_cache))

liouvillian(H::QuantumObject{MT1,SuperOperatorQuantumObject}, Id_cache::Diagonal) where {MT1<:AbstractMatrix} = H

function liouvillian_floquet(
    L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real;
    n_max::Int = 3,
    tol::Real = 1e-15,
) where {T1,T2,T3}
    ((L₀.dims == Lₚ.dims) && (L₀.dims == Lₘ.dims)) ||
        throw(ErrorException("The operators are not of the same Hilbert dimension."))

    return _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, n_max, tol)
end

function liouvillian_floquet(
    H::QuantumObject{<:AbstractArray{T1},OpType1},
    c_ops::AbstractVector,
    Hₚ::QuantumObject{<:AbstractArray{T2},OpType2},
    Hₘ::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real;
    n_max::Int = 3,
    tol::Real = 1e-15,
) where {
    T1,
    T2,
    T3,
    OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
}
    return liouvillian_floquet(liouvillian(H, c_ops), liouvillian(Hₚ), liouvillian(Hₘ), ω, n_max = n_max, tol = tol)
end

@doc raw"""
    liouvillian_generalized(H::QuantumObject, fields::Vector, 
    T_list::Vector; N_trunc::Int=size(H,1), tol::Float64=0.0, σ_filter::Union{Nothing, Real}=nothing)

Constructs the generalized Liouvillian for a system coupled to a bath of harmonic oscillators.

See, e.g., Settineri, Alessio, et al. "Dissipation and thermal noise in hybrid quantum systems in the ultrastrong-coupling regime." Physical Review A 98.5 (2018): 053834.
"""
function liouvillian_generalized(
    H::QuantumObject{MT,OperatorQuantumObject},
    fields::Vector,
    T_list::Vector{<:Real};
    N_trunc::Int = size(H, 1),
    tol::Real = 1e-12,
    σ_filter::Union{Nothing,Real} = nothing,
) where {MT<:AbstractMatrix}
    (length(fields) == length(T_list)) || throw(DimensionMismatch("The number of fields, ωs and Ts must be the same."))

    dims = N_trunc == size(H, 1) ? H.dims : [N_trunc]
    result = eigen(H)
    E = real.(result.values[1:N_trunc])
    U = QuantumObject(result.vectors, result.type, result.dims)

    H_d = QuantumObject(Diagonal(complex(E)), dims = dims)

    Ω = E' .- E
    Ωp = triu(dense_to_sparse(Ω, tol), 1)

    # Filter in the Hilbert space
    σ = isnothing(σ_filter) ? 500 * maximum([norm(field) / length(field) for field in fields]) : σ_filter
    F1 = QuantumObject(gaussian.(Ω, 0, σ), dims = dims)
    F1 = dense_to_sparse(F1, tol)

    # Filter in the Liouville space
    # M1 = ones(N_trunc, N_trunc)
    M1 = similar(E, N_trunc, N_trunc)
    M1 .= 1
    Ω1 = kron(Ω, M1)
    Ω2 = kron(M1, Ω)
    Ωdiff = Ω1 .- Ω2
    F2 = QuantumObject(gaussian.(Ωdiff, 0, σ), SuperOperator, dims)
    F2 = dense_to_sparse(F2, tol)

    L = liouvillian(H_d)

    for i in eachindex(fields)
        # The operator that couples the system to the bath in the eigenbasis
        X_op = dense_to_sparse((U'*fields[i]*U).data[1:N_trunc, 1:N_trunc], tol)
        if ishermitian(fields[i])
            X_op = (X_op + X_op') / 2 # Make sure it's hermitian
        end

        # Ohmic reservoir
        N_th = n_th.(Ωp, T_list[i])
        Sp₀ = QuantumObject(triu(X_op, 1), dims = dims)
        Sp₁ = QuantumObject(droptol!((@. Ωp * N_th * Sp₀.data), tol), dims = dims)
        Sp₂ = QuantumObject(droptol!((@. Ωp * (1 + N_th) * Sp₀.data), tol), dims = dims)
        # S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

        L +=
            1 / 2 *
            (F2 .* (sprepost(Sp₁', Sp₀) + sprepost(Sp₀', Sp₁)) - spre(F1 .* (Sp₀ * Sp₁')) - spost(F1 .* (Sp₁ * Sp₀')))
        L +=
            1 / 2 *
            (F2 .* (sprepost(Sp₂, Sp₀') + sprepost(Sp₀, Sp₂')) - spre(F1 .* (Sp₀' * Sp₂)) - spost(F1 .* (Sp₂' * Sp₀)))
    end

    return E, U, L
end

function _liouvillian_floquet(
    L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real,
    n_max::Int,
    tol::Real,
) where {T1,T2,T3}
    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data
    L_p_dense = sparse_to_dense(Lₚ.data)
    L_m_dense = sparse_to_dense(Lₘ.data)

    S = -(L_0 - 1im * n_max * ω * I) \ L_p_dense
    T = -(L_0 + 1im * n_max * ω * I) \ L_m_dense

    for n_i in n_max-1:-1:1
        S = -(L_0 - 1im * n_i * ω * I + L_m * S) \ L_p_dense
        T = -(L_0 + 1im * n_i * ω * I + L_p * T) \ L_m_dense
    end

    tol == 0 && return QuantumObject(L_0 + L_m * S + L_p * T, SuperOperator, L₀.dims)
    return QuantumObject(dense_to_sparse(L_0 + L_m * S + L_p * T, tol), SuperOperator, L₀.dims)
end
