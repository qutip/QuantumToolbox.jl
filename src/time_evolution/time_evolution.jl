export TimeDependentOperatorSum
export TimeEvolutionSol, TimeEvolutionMCSol

export liouvillian, liouvillian_floquet, liouvillian_generalized

const DEFAULT_ODE_SOLVER_OPTIONS = (abstol = 1e-8, reltol = 1e-6, save_everystep = false, save_end = true)
const DEFAULT_SDE_SOLVER_OPTIONS = (abstol = 1e-2, reltol = 1e-2, save_everystep = false, save_end = true)

@doc raw"""
    struct TimeEvolutionSol

A structure storing the results and some information from solving time evolution.

# Fields (Attributes)

- `times::AbstractVector`: The time list of the evolution.
- `states::Vector{QuantumObject}`: The list of result states.
- `expect::Matrix`: The expectation values corresponding to each time point in `times`.
- `retcode`: The return code from the solver.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionSol{TT<:Vector{<:Real},TS<:AbstractVector,TE<:Matrix{ComplexF64}}
    times::TT
    states::TS
    expect::TE
    retcode::Enum
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm
    abstol::Real
    reltol::Real
end

function Base.show(io::IO, sol::TimeEvolutionSol)
    print(io, "Solution of time evolution\n")
    print(io, "(return code: $(sol.retcode))\n")
    print(io, "--------------------------\n")
    print(io, "num_states = $(length(sol.states))\n")
    print(io, "num_expect = $(size(sol.expect, 1))\n")
    print(io, "ODE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

@doc raw"""
    struct TimeEvolutionMCSol

A structure storing the results and some information from solving quantum trajectories of the Monte Carlo wave function time evolution.

# Fields (Attributes)

- `n_traj::Int`: Number of trajectories
- `times::AbstractVector`: The time list of the evolution in each trajectory.
- `states::Vector{Vector{QuantumObject}}`: The list of result states in each trajectory.
- `expect::Matrix`: The expectation values (averaging all trajectories) corresponding to each time point in `times`.
- `expect_all::Array`: The expectation values corresponding to each trajectory and each time point in `times`
- `jump_times::Vector{Vector{Real}}`: The time records of every quantum jump occurred in each trajectory.
- `jump_which::Vector{Vector{Int}}`: The indices of the jump operators in `c_ops` that describe the corresponding quantum jumps occurred in each trajectory.
- `converged::Bool`: Whether the solution is converged or not.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionMCSol{
    TT<:Vector{<:Vector{<:Real}},
    TS<:AbstractVector,
    TE<:Matrix{ComplexF64},
    TEA<:Array{ComplexF64,3},
    TJT<:Vector{<:Vector{<:Real}},
    TJW<:Vector{<:Vector{<:Integer}},
}
    n_traj::Int
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    jump_times::TJT
    jump_which::TJW
    converged::Bool
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm
    abstol::Real
    reltol::Real
end

function Base.show(io::IO, sol::TimeEvolutionMCSol)
    print(io, "Solution of quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.n_traj)\n")
    print(io, "num_states = $(length(sol.states[1]))\n")
    print(io, "num_expect = $(size(sol.expect, 1))\n")
    print(io, "ODE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

struct TimeEvolutionSSESol{
    TT<:Vector{<:Real},
    TS<:AbstractVector,
    TE<:Matrix{ComplexF64},
    TEA<:Array{ComplexF64,3},
    T1<:Real,
    T2<:Real,
}
    n_traj::Int
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    converged::Bool
    alg::StochasticDiffEq.StochasticDiffEqAlgorithm
    abstol::T1
    reltol::T2
end

function Base.show(io::IO, sol::TimeEvolutionSSESol)
    print(io, "Solution of quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.n_traj)\n")
    # print(io, "num_states = $(length(sol.states[1]))\n")
    print(io, "num_expect = $(size(sol.expect, 1))\n")
    print(io, "SDE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

abstract type LindbladJumpCallbackType end

struct ContinuousLindbladJumpCallback <: LindbladJumpCallbackType
    interp_points::Int
end

struct DiscreteLindbladJumpCallback <: LindbladJumpCallbackType end

ContinuousLindbladJumpCallback(; interp_points::Int = 10) = ContinuousLindbladJumpCallback(interp_points)

## Time-dependent sum of operators

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
    liouvillian(H::QuantumObject, c_ops::AbstractVector, Id_cache=I(prod(H.dims)))

Construct the Liouvillian [`SuperOperator`](@ref) for a system Hamiltonian ``\hat{H}`` and a set of collapse operators ``\{\hat{C}_n\}_n``:

```math
\mathcal{L} [\cdot] = -i[\hat{H}, \cdot] + \sum_n \mathcal{D}(\hat{C}_n) [\cdot]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\cdot] = \hat{C}_n [\cdot] \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n [\cdot] - \frac{1}{2} [\cdot] \hat{C}_n^\dagger \hat{C}_n
```

The optional argument `Id_cache` can be used to pass a precomputed identity matrix. This can be useful when the same function is applied multiple times with a known Hilbert space dimension.

See also [`spre`](@ref), [`spost`](@ref), and [`lindblad_dissipator`](@ref).
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
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))

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
