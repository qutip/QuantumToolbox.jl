export TimeEvolutionSol
export TimeEvolutionMultiTrajSol, TimeEvolutionMCSol, TimeEvolutionStochasticSol
export average_states, average_expect, std_expect

export liouvillian_floquet, liouvillian_dressed_nonsecular

default_ode_solver_options(::Type{T}) where {T} = (abstol = _float_type(T)(1.0e-8), reltol = _float_type(T)(1.0e-6), save_everystep = false, save_end = true)
default_sde_solver_options(::Type{T}) where {T} = (abstol = _float_type(T)(1.0e-3), reltol = _float_type(T)(2.0e-3), save_everystep = false, save_end = true)
const COL_TIMES_WHICH_INIT_SIZE = 200

abstract type TimeEvolutionMultiTrajSol{Tstates, Texpect} end

@doc raw"""
    struct TimeEvolutionProblem

A Julia constructor for handling the `ODEProblem` of the time evolution of quantum systems.

# Fields (Attributes)

- `prob::AbstractSciMLProblem`: The `ODEProblem` of the time evolution.
- `times::AbstractVector`: The time list of the evolution.
- `states_type::QuantumObjectType`: The type of the quantum states during the evolution (e.g., [`Ket`](@ref), [`Operator`](@ref), [`OperatorKet`](@ref), or [`SuperOperator`](@ref)).
- `dimensions::AbstractDimensions`: The dimensions of the Hilbert space.
- `kwargs::KWT`: Generic keyword arguments.

!!! note "`dims` property"
    For a given `prob::TimeEvolutionProblem`, `prob.dims` or `getproperty(prob, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct TimeEvolutionProblem{
        ST <: QuantumObjectType,
        DT <: AbstractDimensions,
        PT <: AbstractSciMLProblem,
        TT <: AbstractVector,
        KWT,
    }
    prob::PT
    times::TT
    states_type::ST
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

TimeEvolutionProblem(prob, times, states_type, dims) = TimeEvolutionProblem(prob, times, states_type, dims, nothing)

@doc raw"""
    struct TimeEvolutionSol

A structure storing the results and some information from solving time evolution.

# Fields (Attributes)

- `times::AbstractVector`: The list of time points at which the expectation values are calculated during the evolution.
- `times_states::AbstractVector`: The list of time points at which the states are stored during the evolution.
- `states::Vector{QuantumObject}`: The list of result states corresponding to each time point in `times_states`.
- `expect::Union{AbstractMatrix,Nothing}`: The expectation values corresponding to each time point in `times`.
- `retcode`: The return code from the solver.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionSol{
        TT1 <: AbstractVector{<:Real},
        TT2 <: AbstractVector{<:Real},
        TS <: AbstractVector,
        TE <: Union{AbstractMatrix, Nothing},
        RETT <: Enum,
        AlgT <: AbstractODEAlgorithm,
        TolT <: Real,
    }
    times::TT1
    times_states::TT2
    states::TS
    expect::TE
    retcode::RETT
    alg::AlgT
    abstol::TolT
    reltol::TolT
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
    struct TimeEvolutionMCSol <: TimeEvolutionMultiTrajSol

A structure storing the results and some information from solving quantum trajectories of the Monte Carlo wave function time evolution.

# Fields (Attributes)

- `ntraj::Int`: Number of trajectories
- `times::AbstractVector`: The list of time points at which the expectation values are calculated during the evolution.
- `times_states::AbstractVector`: The list of time points at which the states are stored during the evolution.
- `states::AbstractVecOrMat`: The list of result states in each trajectory and each time point in `times_states`.
- `expect::Union{AbstractArray,Nothing}`: The expectation values corresponding to each trajectory and each time point in `times`.
- `col_times::Vector{Vector{Real}}`: The time records of every quantum jump occurred in each trajectory.
- `col_which::Vector{Vector{Int}}`: The indices of which collapse operator was responsible for each quantum jump in `col_times`.
- `converged::Bool`: Whether the solution is converged or not.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.

# Notes

When the keyword argument `keep_runs_results` is passed as `Val(false)` to a multi-trajectory solver, the `states` and `expect` fields store only the average results (averaged over all trajectories). The results can be accessed by the following index-order:

- `sol.states[time_idx]`
- `sol.expect[e_op,time_idx]`

If the keyword argument `keep_runs_results = Val(true)`, the results for each trajectory and each time are stored, and the index-order of the elements in fields `states` and `expect` are:

- `sol.states[trajectory,time_idx]`
- `sol.expect[e_op,trajectory,time_idx]`

We also provide the following functions for statistical analysis of multi-trajectory solutions:

- [`average_states`](@ref)
- [`average_expect`](@ref)
- [`std_expect`](@ref)
"""
struct TimeEvolutionMCSol{
        TT1 <: AbstractVector{<:Real},
        TT2 <: AbstractVector{<:Real},
        TS <: AbstractVecOrMat,
        TE <: Union{AbstractArray, Nothing},
        TJT <: Vector{<:Vector{<:Real}},
        TJW <: Vector{<:Vector{<:Integer}},
        AlgT <: AbstractODEAlgorithm,
        TolT <: Real,
    } <: TimeEvolutionMultiTrajSol{TS, TE}
    ntraj::Int
    times::TT1
    times_states::TT2
    states::TS
    expect::TE
    col_times::TJT
    col_which::TJW
    converged::Bool
    alg::AlgT
    abstol::TolT
    reltol::TolT
end

function Base.show(io::IO, sol::TimeEvolutionMCSol)
    print(io, "Solution of quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.ntraj)\n")
    print(io, "num_states = $(size(sol.states, ndims(sol.states)))\n") # get the size of last dimension
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
    struct TimeEvolutionStochasticSol

A structure storing the results and some information from solving trajectories of the Stochastic time evolution.

# Fields (Attributes)

- `ntraj::Int`: Number of trajectories
- `times::AbstractVector`: The list of time points at which the expectation values are calculated during the evolution.
- `times_states::AbstractVector`: The list of time points at which the states are stored during the evolution.
- `states::AbstractVecOrMat`: The list of result states in each trajectory and each time point in `times_states`.
- `expect::Union{AbstractArray,Nothing}`: The expectation values corresponding to each trajectory and each time point in `times`.
- `converged::Bool`: Whether the solution is converged or not.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.

# Notes

When the keyword argument `keep_runs_results` is passed as `Val(false)` to a multi-trajectory solver, the `states` and `expect` fields store only the average results (averaged over all trajectories). The results can be accessed by the following index-order:

- `sol.states[time_idx]`
- `sol.expect[e_op,time_idx]`

If the keyword argument `keep_runs_results = Val(true)`, the results for each trajectory and each time are stored, and the index-order of the elements in fields `states` and `expect` are:

- `sol.states[trajectory,time_idx]`
- `sol.expect[e_op,trajectory,time_idx]`

We also provide the following functions for statistical analysis of multi-trajectory solutions:

- [`average_states`](@ref)
- [`average_expect`](@ref)
- [`std_expect`](@ref)
"""
struct TimeEvolutionStochasticSol{
        TT1 <: AbstractVector{<:Real},
        TT2 <: AbstractVector{<:Real},
        TS <: AbstractVecOrMat,
        TE <: Union{AbstractArray, Nothing},
        TEM <: Union{AbstractArray, Nothing},
        AlgT <: AbstractSDEAlgorithm,
        TolT <: Real,
    } <: TimeEvolutionMultiTrajSol{TS, TE}
    ntraj::Int
    times::TT1
    times_states::TT2
    states::TS
    expect::TE
    measurement::TEM
    converged::Bool
    alg::AlgT
    abstol::TolT
    reltol::TolT
end

function Base.show(io::IO, sol::TimeEvolutionStochasticSol)
    print(io, "Solution of stochastic quantum trajectories\n")
    print(io, "(converged: $(sol.converged))\n")
    print(io, "--------------------------------\n")
    print(io, "num_trajectories = $(sol.ntraj)\n")
    print(io, "num_states = $(size(sol.states, ndims(sol.states)))\n") # get the size of last dimension
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

@doc raw"""
    average_states(sol::TimeEvolutionMultiTrajSol)

Return the trajectory-averaged result states (as density [`Operator`](@ref)) at each time point.
"""
average_states(sol::TimeEvolutionMultiTrajSol{<:Matrix{<:QuantumObject}}) = _average_traj_states(sol.states)
average_states(sol::TimeEvolutionMultiTrajSol{<:Vector{<:QuantumObject}}) = sol.states  # this case should already be averaged over all trajectories

# TODO: Check if broadcasting division ./ size(states, 1) is type stable
_average_traj_states(states::Matrix{<:QuantumObject{Ket}}) =
    map(x -> x / size(states, 1), dropdims(sum(ket2dm, states, dims = 1), dims = 1))
_average_traj_states(states::Matrix{<:QuantumObject{ObjType}}) where {ObjType <: Union{Operator, OperatorKet}} =
    map(x -> x / size(states, 1), dropdims(sum(states, dims = 1), dims = 1))

@doc raw"""
    average_expect(sol::TimeEvolutionMultiTrajSol)

Return the trajectory-averaged expectation values at each time point.
"""
average_expect(sol::TimeEvolutionMultiTrajSol{TS, Array{T, 3}}) where {TS, T <: Number} = _average_traj_expect(sol.expect)
average_expect(sol::TimeEvolutionMultiTrajSol{TS, Matrix{T}}) where {TS, T <: Number} = sol.expect  # this case should already be averaged over all trajectories
average_expect(::TimeEvolutionMultiTrajSol{TS, Nothing}) where {TS} = nothing

_average_traj_expect(expvals::Array{T, 3}) where {T <: Number} =
    dropdims(sum(expvals, dims = 2), dims = 2) ./ size(expvals, 2)

# these are used in multi-trajectory solvers before returning solutions
_store_multitraj_states(states::Matrix{<:QuantumObject}, keep_runs_results::Val{false}) = _average_traj_states(states)
_store_multitraj_states(states::Matrix{<:QuantumObject}, keep_runs_results::Val{true}) = states
_store_multitraj_expect(expvals::Array{T, 3}, keep_runs_results::Val{false}) where {T <: Number} =
    _average_traj_expect(expvals)
_store_multitraj_expect(expvals::Array{T, 3}, keep_runs_results::Val{true}) where {T <: Number} = expvals
_store_multitraj_expect(expvals::Nothing, keep_runs_results) = nothing

@doc raw"""
    std_expect(sol::TimeEvolutionMultiTrajSol)

Return the trajectory-wise standard deviation of the expectation values at each time point.
"""
function std_expect(sol::TimeEvolutionMultiTrajSol{TS, Array{T, 3}}) where {TS, T <: Number}
    # the following standard deviation (std) is defined as the square-root of variance instead of pseudo-variance
    # i.e., it is equivalent to (even for complex expectation values):
    #    dropdims(
    #        sqrt.(mean(abs2.(sol.expect), dims = 2) .- abs2.(mean(sol.expect, dims = 2))),
    #        dims = 2
    #    )
    # [this should be included in the runtest]
    return dropdims(std(sol.expect, corrected = false, dims = 2), dims = 2)
end
std_expect(::TimeEvolutionMultiTrajSol{TS, Matrix{T}}) where {TS, T <: Number} = throw(
    ArgumentError(
        "Can not compute the standard deviation without the expectation values of each trajectory. Try to specify keyword argument `keep_runs_results=Val(true)` to the solver.",
    ),
)
std_expect(::TimeEvolutionMultiTrajSol{TS, Nothing}) where {TS} = nothing

#######################################
#=
    Callbacks for Monte Carlo quantum trajectories
=#

abstract type LindbladJumpCallbackType end

struct ContinuousLindbladJumpCallback <: LindbladJumpCallbackType
    interp_points::Int
end

struct DiscreteLindbladJumpCallback <: LindbladJumpCallbackType end

ContinuousLindbladJumpCallback(; interp_points::Int = 10) = ContinuousLindbladJumpCallback(interp_points)

function _check_tlist(tlist, T::Type)
    tlist2 = convert(Vector{T}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    # Check if the list of times is not empty
    isempty(tlist2) && throw(ArgumentError("The time list must not be empty."))
    # Check if the list of times is sorted
    issorted(tlist2) || throw(ArgumentError("The time list must be sorted."))
    # Check if the list of times is unique
    allunique(tlist2) || throw(ArgumentError("The time list must be unique."))

    return tlist2
end

#######################################

_make_c_ops_list(c_ops) = c_ops
_make_c_ops_list(c_ops::AbstractQuantumObject) = (c_ops,)

function _merge_saveat(tlist, e_ops, default_options; kwargs...)
    is_empty_e_ops = isnothing(e_ops) ? true : isempty(e_ops)
    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (default_options..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)

    # DifferentialEquations.jl has this weird save_end setting
    # So we need to do this to make sure it's consistent
    haskey(kwargs, :save_end) && return kwargs2
    isempty(kwargs2.saveat) && return kwargs2

    save_end = tlist[end] in kwargs2.saveat
    return merge(kwargs2, (save_end = save_end,))
end

#######################################
#=
Helpers for handling output of ensemble problems.
This is very useful especially for dispatching which method to use to update the progress bar.
=#

# Output function with progress bar update
function _ensemble_output_func_progress(sol, i, progr, output_func)
    next!(progr)
    return output_func(sol, i)
end

# Output function with distributed channel update for progress bar
function _ensemble_output_func_distributed(sol, i, channel, output_func)
    put!(channel, true)
    return output_func(sol, i)
end

function _ensemble_dispatch_output_func(
        ::ET,
        progress_bar,
        ntraj,
        output_func;
        progr_desc = "Progress: ",
    ) where {ET <: Union{EnsembleSerial, EnsembleThreads}}
    if getVal(progress_bar)
        progr = Progress(ntraj; enabled = getVal(progress_bar), desc = progr_desc, settings.ProgressMeterKWARGS...)
        f = (sol, i) -> _ensemble_output_func_progress(sol, i, progr, output_func)
        return (f, progr, nothing)
    else
        return (output_func, nothing, nothing)
    end
end
function _ensemble_dispatch_output_func(
        ::ET,
        progress_bar,
        ntraj,
        output_func;
        progr_desc = "Progress... ",
    ) where {ET <: Union{EnsembleSplitThreads, EnsembleDistributed}}
    if getVal(progress_bar)
        progr = Progress(ntraj; enabled = getVal(progress_bar), desc = progr_desc, settings.ProgressMeterKWARGS...)
        progr_channel::RemoteChannel{Channel{Bool}} = RemoteChannel(() -> Channel{Bool}(1))

        f = (sol, i) -> _ensemble_output_func_distributed(sol, i, progr_channel, output_func)
        return (f, progr, progr_channel)
    else
        return (output_func, nothing, nothing)
    end
end

function _ensemble_dispatch_prob_func(rng, ntraj, tlist, prob_func; kwargs...)
    seeds = map(i -> rand(rng, UInt64), 1:ntraj)
    return (prob, i, repeat) -> prob_func(prob, i, repeat, rng, seeds, tlist; kwargs...)
end

function _ensemble_dispatch_solve(
        ens_prob::TimeEvolutionProblem,
        alg::Union{<:AbstractODEAlgorithm, <:AbstractSDEAlgorithm},
        ensemblealg::ET,
        ntraj::Int,
    ) where {ET <: Union{EnsembleSplitThreads, EnsembleDistributed}}
    sol = nothing

    @sync begin
        @async while take!(ens_prob.kwargs.channel)
            next!(ens_prob.kwargs.progr)
        end

        @async begin
            sol = solve(ens_prob.prob, alg, ensemblealg, trajectories = ntraj)
            put!(ens_prob.kwargs.channel, false)
        end
    end

    return sol
end
function _ensemble_dispatch_solve(
        ens_prob::TimeEvolutionProblem,
        alg::Union{<:AbstractODEAlgorithm, <:AbstractSDEAlgorithm},
        ensemblealg,
        ntraj::Int,
    )
    sol = solve(ens_prob.prob, alg, ensemblealg, trajectories = ntraj)
    return sol
end

# For mapped solvers
function _se_me_map_prob_func(prob, i, repeat, iter)
    f = deepcopy(prob.f.f)
    u0 = iter[i][1]
    p = iter[i][2:end]
    if haskey(prob.kwargs, :callback)
        return remake(prob, f = f, u0 = u0, p = p, callback = deepcopy(prob.kwargs[:callback]))
    else
        return remake(prob, f = f, u0 = u0, p = p)
    end
end

#######################################
#=
 Stochastic funcs
=#
function _stochastic_prob_func(prob, i, repeat, rng, seeds, tlist; kwargs...)
    seed = seeds[i]
    traj_rng = typeof(rng)()
    seed!(traj_rng, seed)

    sc_ops = kwargs[:sc_ops]
    store_measurement = kwargs[:store_measurement]
    noise = _make_noise(prob.prob.tspan[1], sc_ops, store_measurement, traj_rng)

    return remake(prob.prob, noise = noise, seed = seed)
end

# Standard output function which does nothing (used for mapped and stochastic solvers)
_standard_output_func(sol, i) = (sol, false)

#= 
    Define diagonal or non-diagonal noise depending on the type of `sc_ops`.
    If `sc_ops` is a `AbstractQuantumObject`, we avoid using the non-diagonal noise.
=#
function _make_noise(t0, sc_ops, store_measurement::Val, rng)
    noise = RealWienerProcess!(
        t0,
        zeros(length(sc_ops)),
        zeros(length(sc_ops)),
        save_everystep = getVal(store_measurement),
        rng = rng,
    )

    return noise
end
function _make_noise(t0, sc_ops::AbstractQuantumObject, store_measurement::Val, rng)
    noise = RealWienerProcess(t0, 0.0, 0.0, save_everystep = getVal(store_measurement), rng = rng)

    return noise
end

#=
    struct DiffusionOperator

A struct to represent the diffusion operator. This is used to perform the diffusion process on N different Wiener processes.
=#
struct DiffusionOperator{T, OpType <: Union{Tuple{Vararg{AbstractSciMLOperator}}, AbstractVector{<:AbstractSciMLOperator}}}
    ops::OpType
    function DiffusionOperator(ops::OpType) where {OpType}
        T = mapreduce(eltype, promote_type, ops)
        return new{T, OpType}(ops)
    end
end

@generated function (L::DiffusionOperator{T, OpType} where {T, OpType <: Tuple})(w, v, p, t)
    ops_types = L.parameters[2].parameters
    N = length(ops_types)
    return quote
        M = length(v)
        S = (size(w, 1), size(w, 2)) # This supports also `w` as a `Vector`
        (S[1] == M && S[2] == $N) || throw(DimensionMismatch("The size of the output vector is incorrect."))
        Base.@nexprs $N i -> begin
            L.ops[i](@view(w[:, i]), v, v, p, t)
        end
        return w
    end
end

function (L::DiffusionOperator{T, OpType} where {T, OpType <: AbstractVector})(w, v, p, t)
    N = length(L.ops)
    M = length(v)
    S = (size(w, 1), size(w, 2)) # This supports also `w` as a `Vector`
    (S[1] == M && S[2] == N) || throw(DimensionMismatch("The size of the output vector is incorrect."))

    for i in Base.OneTo(N)
        L.ops[i](@view(w[:, i]), v, v, p, t)
    end
    return w
end

#######################################

function liouvillian_floquet(
        L₀::QuantumObject{SuperOperator},
        Lₚ::QuantumObject{SuperOperator},
        Lₘ::QuantumObject{SuperOperator},
        ω::Real;
        n_max::Int = 3,
        tol::Real = 1.0e-15,
    )
    check_dimensions(L₀, Lₚ, Lₘ)
    return _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, n_max, tol)
end

function liouvillian_floquet(
        H::QuantumObject{OpType1},
        Hₚ::QuantumObject{OpType2},
        Hₘ::QuantumObject{OpType3},
        ω::Real,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        n_max::Int = 3,
        tol::Real = 1.0e-15,
    ) where {
        OpType1 <: Union{Operator, SuperOperator},
        OpType2 <: Union{Operator, SuperOperator},
        OpType3 <: Union{Operator, SuperOperator},
    }
    return liouvillian_floquet(liouvillian(H, c_ops), liouvillian(Hₚ), liouvillian(Hₘ), ω, n_max = n_max, tol = tol)
end

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
    (size(A, 1) == N && size(B, 1) == N) || throw(DimensionMismatch("Matrix sizes do not match energy list."))

    I_A, J_A, V_A = _findnz(A)
    I_B, J_B, V_B = _findnz(B)

    nnzA = length(V_A)
    nnzB = length(V_B)

    T = Base.promote_eltype(V_A, V_B, E, σ)
    I_out = Int[]
    J_out = Int[]
    V_out = Vector{T}()

    nnzA == 0 || nnzB == 0 && return sparse(I_out, J_out, V_out, N^2, N^2)

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

function _liouvillian_floquet(
        L₀::QuantumObject{SuperOperator},
        Lₚ::QuantumObject{SuperOperator},
        Lₘ::QuantumObject{SuperOperator},
        ω::Real,
        n_max::Int,
        tol::Real,
    )
    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data
    L_p_dense = to_dense(Lₚ.data)
    L_m_dense = to_dense(Lₘ.data)

    S = -(L_0 - 1im * n_max * ω * I) \ L_p_dense
    T = -(L_0 + 1im * n_max * ω * I) \ L_m_dense

    for n_i in (n_max - 1):-1:1
        S = -(L_0 - 1im * n_i * ω * I + L_m * S) \ L_p_dense
        T = -(L_0 + 1im * n_i * ω * I + L_p * T) \ L_m_dense
    end

    tol == 0 && return QuantumObject(L_0 + L_m * S + L_p * T, SuperOperator(), L₀.dimensions)
    return QuantumObject(to_sparse(L_0 + L_m * S + L_p * T, tol), SuperOperator(), L₀.dimensions)
end
