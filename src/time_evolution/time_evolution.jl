export TimeEvolutionSol, TimeEvolutionMCSol, TimeEvolutionSSESol

export liouvillian_floquet, liouvillian_generalized

const DEFAULT_ODE_SOLVER_OPTIONS = (abstol = 1e-8, reltol = 1e-6, save_everystep = false, save_end = true)
const DEFAULT_SDE_SOLVER_OPTIONS = (abstol = 1e-2, reltol = 1e-2, save_everystep = false, save_end = true)
const JUMP_TIMES_WHICH_INIT_SIZE = 200

@doc raw"""
    struct TimeEvolutionProblem

A Julia constructor for handling the `ODEProblem` of the time evolution of quantum systems.

# Fields (Attributes)

- `prob::AbstractSciMLProblem`: The `ODEProblem` of the time evolution.
- `times::Abstractvector`: The time list of the evolution.
- `dimensions::AbstractDimensions`: The dimensions of the Hilbert space.
- `kwargs::KWT`: Generic keyword arguments.

!!! note "`dims` property"
    For a given `prob::TimeEvolutionProblem`, `prob.dims` or `getproperty(prob, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct TimeEvolutionProblem{PT<:AbstractSciMLProblem,TT<:AbstractVector,DT<:AbstractDimensions,KWT}
    prob::PT
    times::TT
    dimensions::DT
    kwargs::KWT
end

function Base.getproperty(prob::TimeEvolutionProblem, key::Symbol)
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(prob, :dimensions))
    else
        return getfield(prob, key)
    end
end

TimeEvolutionProblem(prob, times, dims) = TimeEvolutionProblem(prob, times, dims, nothing)

@doc raw"""
    struct TimeEvolutionSol

A structure storing the results and some information from solving time evolution.

# Fields (Attributes)

- `times::AbstractVector`: The time list of the evolution.
- `states::Vector{QuantumObject}`: The list of result states.
- `expect::Union{AbstractMatrix,Nothing}`: The expectation values corresponding to each time point in `times`.
- `retcode`: The return code from the solver.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionSol{
    TT<:AbstractVector{<:Real},
    TS<:AbstractVector,
    TE<:Union{AbstractMatrix,Nothing},
    RETT<:Enum,
    AlgT<:OrdinaryDiffEqAlgorithm,
    AT<:Real,
    RT<:Real,
}
    times::TT
    states::TS
    expect::TE
    retcode::RETT
    alg::AlgT
    abstol::AT
    reltol::RT
end

function Base.show(io::IO, sol::TimeEvolutionSol)
    print(io, "Solution of time evolution\n")
    print(io, "(return code: $(sol.retcode))\n")
    print(io, "--------------------------\n")
    print(io, "num_states = $(length(sol.states))\n")
    if sol.expect isa Nothing
        print(io, "num_expect = 0\n")
    else
        print(io, "num_expect = $(size(sol.expect, 1))\n")
    end
    print(io, "ODE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

@doc raw"""
    struct TimeEvolutionMCSol

A structure storing the results and some information from solving quantum trajectories of the Monte Carlo wave function time evolution.

# Fields (Attributes)

- `ntraj::Int`: Number of trajectories
- `times::AbstractVector`: The time list of the evolution.
- `states::Vector{Vector{QuantumObject}}`: The list of result states in each trajectory.
- `expect::Union{AbstractMatrix,Nothing}`: The expectation values (averaging all trajectories) corresponding to each time point in `times`.
- `expect_all::Union{AbstractMatrix,Nothing}`: The expectation values corresponding to each trajectory and each time point in `times`
- `jump_times::Vector{Vector{Real}}`: The time records of every quantum jump occurred in each trajectory.
- `jump_which::Vector{Vector{Int}}`: The indices of the jump operators in `c_ops` that describe the corresponding quantum jumps occurred in each trajectory.
- `converged::Bool`: Whether the solution is converged or not.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionMCSol{
    TT<:AbstractVector{<:Real},
    TS<:AbstractVector,
    TE<:Union{AbstractMatrix,Nothing},
    TEA<:Union{AbstractArray,Nothing},
    TJT<:Vector{<:Vector{<:Real}},
    TJW<:Vector{<:Vector{<:Integer}},
    AlgT<:OrdinaryDiffEqAlgorithm,
    AT<:Real,
    RT<:Real,
}
    ntraj::Int
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    jump_times::TJT
    jump_which::TJW
    converged::Bool
    alg::AlgT
    abstol::AT
    reltol::RT
end

function Base.show(io::IO, sol::TimeEvolutionMCSol)
    print(io, "Solution of quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.ntraj)\n")
    print(io, "num_states = $(length(sol.states[1]))\n")
    if sol.expect isa Nothing
        print(io, "num_expect = 0\n")
    else
        print(io, "num_expect = $(size(sol.expect, 1))\n")
    end
    print(io, "ODE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

@doc raw"""
    struct TimeEvolutionSSESol

A structure storing the results and some information from solving trajectories of the Stochastic Shrodinger equation time evolution.

# Fields (Attributes)

- `ntraj::Int`: Number of trajectories
- `times::AbstractVector`: The time list of the evolution.
- `states::Vector{Vector{QuantumObject}}`: The list of result states in each trajectory.
- `expect::Union{AbstractMatrix,Nothing}`: The expectation values (averaging all trajectories) corresponding to each time point in `times`.
- `expect_all::Union{AbstractArray,Nothing}`: The expectation values corresponding to each trajectory and each time point in `times`
- `converged::Bool`: Whether the solution is converged or not.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionSSESol{
    TT<:AbstractVector{<:Real},
    TS<:AbstractVector,
    TE<:Union{AbstractMatrix,Nothing},
    TEA<:Union{AbstractArray,Nothing},
    AlgT<:StochasticDiffEqAlgorithm,
    AT<:Real,
    RT<:Real,
}
    ntraj::Int
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    converged::Bool
    alg::AlgT
    abstol::AT
    reltol::RT
end

function Base.show(io::IO, sol::TimeEvolutionSSESol)
    print(io, "Solution of quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.ntraj)\n")
    print(io, "num_states = $(length(sol.states[1]))\n")
    if sol.expect isa Nothing
        print(io, "num_expect = 0\n")
    else
        print(io, "num_expect = $(size(sol.expect, 1))\n")
    end
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

function _check_tlist(tlist, T::Type)
    tlist2 = convert(Vector{T}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    # Check if the list of times is not empty
    isempty(tlist2) && throw(ArgumentError("The list of times must not be empty."))
    # Check if the list of times is sorted
    !issorted(tlist2) && throw(ArgumentError("The list of times must be sorted."))
    # Check if the list of times is unique
    length(tlist2) != length(unique(tlist2)) && throw(ArgumentError("The list of times must be unique."))

    return tlist2
end

#######################################

function liouvillian_floquet(
    L₀::QuantumObject{SuperOperatorQuantumObject},
    Lₚ::QuantumObject{SuperOperatorQuantumObject},
    Lₘ::QuantumObject{SuperOperatorQuantumObject},
    ω::Real;
    n_max::Int = 3,
    tol::Real = 1e-15,
)
    check_dimensions(L₀, Lₚ, Lₘ)
    return _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, n_max, tol)
end

function liouvillian_floquet(
    H::QuantumObject{OpType1},
    Hₚ::QuantumObject{OpType2},
    Hₘ::QuantumObject{OpType3},
    ω::Real,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    n_max::Int = 3,
    tol::Real = 1e-15,
) where {
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
    H::QuantumObject{OperatorQuantumObject},
    fields::Vector,
    T_list::Vector{<:Real};
    N_trunc::Union{Int,Nothing} = nothing,
    tol::Real = 1e-12,
    σ_filter::Union{Nothing,Real} = nothing,
)
    (length(fields) == length(T_list)) || throw(DimensionMismatch("The number of fields, ωs and Ts must be the same."))

    dims = (N_trunc isa Nothing) ? H.dimensions : SVector(N_trunc)
    final_size = prod(dims)
    result = eigen(H)
    E = real.(result.values[1:final_size])
    U = QuantumObject(result.vectors, result.type, result.dimensions)

    H_d = QuantumObject(Diagonal(complex(E)), type = Operator, dims = dims)

    Ω = E' .- E
    Ωp = triu(dense_to_sparse(Ω, tol), 1)

    # Filter in the Hilbert space
    σ = isnothing(σ_filter) ? 500 * maximum([norm(field) / length(field) for field in fields]) : σ_filter
    F1 = QuantumObject(gaussian.(Ω, 0, σ), type = Operator, dims = dims)
    F1 = dense_to_sparse(F1, tol)

    # Filter in the Liouville space
    # M1 = ones(final_size, final_size)
    M1 = similar(E, final_size, final_size)
    M1 .= 1
    Ω1 = kron(Ω, M1)
    Ω2 = kron(M1, Ω)
    Ωdiff = Ω1 .- Ω2
    F2 = QuantumObject(gaussian.(Ωdiff, 0, σ), SuperOperator, dims)
    F2 = dense_to_sparse(F2, tol)

    L = liouvillian(H_d)

    for i in eachindex(fields)
        # The operator that couples the system to the bath in the eigenbasis
        X_op = dense_to_sparse((U'*fields[i]*U).data[1:final_size, 1:final_size], tol)
        if ishermitian(fields[i])
            X_op = (X_op + X_op') / 2 # Make sure it's hermitian
        end

        # Ohmic reservoir
        N_th = n_thermal.(Ωp, T_list[i])
        Sp₀ = QuantumObject(triu(X_op, 1), type = Operator, dims = dims)
        Sp₁ = QuantumObject(droptol!((@. Ωp * N_th * Sp₀.data), tol), type = Operator, dims = dims)
        Sp₂ = QuantumObject(droptol!((@. Ωp * (1 + N_th) * Sp₀.data), tol), type = Operator, dims = dims)
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
    L₀::QuantumObject{SuperOperatorQuantumObject},
    Lₚ::QuantumObject{SuperOperatorQuantumObject},
    Lₘ::QuantumObject{SuperOperatorQuantumObject},
    ω::Real,
    n_max::Int,
    tol::Real,
)
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

    tol == 0 && return QuantumObject(L_0 + L_m * S + L_p * T, SuperOperator, L₀.dimensions)
    return QuantumObject(dense_to_sparse(L_0 + L_m * S + L_p * T, tol), SuperOperator, L₀.dimensions)
end
