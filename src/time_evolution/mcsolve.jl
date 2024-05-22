export mcsolveProblem, mcsolveEnsembleProblem, mcsolve
export ContinuousLindbladJumpCallback, DiscreteLindbladJumpCallback

function _save_func_mcsolve(integrator)
    internal_params = integrator.p
    progr = internal_params.progr_mc

    if !internal_params.is_empty_e_ops_mc
        e_ops = internal_params.e_ops_mc
        expvals = internal_params.expvals
        cache_mc = internal_params.cache_mc

        copyto!(cache_mc, integrator.u)
        normalize!(cache_mc)
        ψ = cache_mc
        _expect = op -> dot(ψ, op, ψ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    next!(progr)
    return u_modified!(integrator, false)
end

function LindbladJumpAffect!(integrator)
    internal_params = integrator.p
    c_ops = internal_params.c_ops
    cache_mc = internal_params.cache_mc
    weights_mc = internal_params.weights_mc
    cumsum_weights_mc = internal_params.cumsum_weights_mc
    random_n = internal_params.random_n
    jump_times = internal_params.jump_times
    jump_which = internal_params.jump_which
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        mul!(cache_mc, c_ops[i], ψ)
        weights_mc[i] = real(dot(cache_mc, cache_mc))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    collaps_idx = getindex(1:length(weights_mc), findfirst(>(rand() * sum(weights_mc)), cumsum_weights_mc))
    mul!(cache_mc, c_ops[collaps_idx], ψ)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)

    # push!(jump_times, integrator.t)
    # push!(jump_which, collaps_idx)
    random_n[] = rand()
    jump_times[internal_params.jump_times_which_idx[]] = integrator.t
    jump_which[internal_params.jump_times_which_idx[]] = collaps_idx
    internal_params.jump_times_which_idx[] += 1
    if internal_params.jump_times_which_idx[] > length(jump_times)
        resize!(jump_times, length(jump_times) + internal_params.jump_times_which_init_size)
        resize!(jump_which, length(jump_which) + internal_params.jump_times_which_init_size)
    end
end

LindbladJumpContinuousCondition(u, t, integrator) = integrator.p.random_n[] - real(dot(u, u))

LindbladJumpDiscreteCondition(u, t, integrator) = real(dot(u, u)) < integrator.p.random_n[]

function _mcsolve_prob_func(prob, i, repeat)
    internal_params = prob.p
    seeds = internal_params.seeds
    !isnothing(seeds) && Random.seed!(seeds[i])

    prm = merge(
        internal_params,
        (
            expvals = similar(internal_params.expvals),
            cache_mc = similar(internal_params.cache_mc),
            weights_mc = similar(internal_params.weights_mc),
            cumsum_weights_mc = similar(internal_params.weights_mc),
            random_n = Ref(rand()),
            progr_mc = ProgressBar(size(internal_params.expvals, 2), enable = false),
            jump_times_which_idx = Ref(1),
            jump_times = similar(internal_params.jump_times),
            jump_which = similar(internal_params.jump_which),
        ),
    )

    return remake(prob, p = prm)
end

function _mcsolve_output_func(sol, i)
    resize!(sol.prob.p.jump_times, sol.prob.p.jump_times_which_idx[] - 1)
    resize!(sol.prob.p.jump_which, sol.prob.p.jump_times_which_idx[] - 1)
    return (sol, false)
end

function _mcsolve_generate_statistics(sol, i, times, states, expvals_all, jump_times, jump_which)
    sol_i = sol[:, i]
    sol_u =
        haskey(sol_i.prob.kwargs, :save_idxs) ? sol_i.u :
        [QuantumObject(sol_i.u[i], dims = sol_i.prob.p.Hdims) for i in 1:length(sol_i.u)]

    copyto!(view(expvals_all, i, :, :), sol_i.prob.p.expvals)
    times[i] = sol_i.t
    states[i] = sol_u
    jump_times[i] = sol_i.prob.p.jump_times
    return jump_which[i] = sol_i.prob.p.jump_which
end

"""
    mcsolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        t_l::AbstractVector,
        c_ops::Vector{QuantumObject{Tc, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[];
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Vector{QuantumObject{Te, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        jump_callback::TJC=ContinuousLindbladJumpCallback(),
        kwargs...)

Generates the ODEProblem for a single trajectory of the Monte Carlo wave function
time evolution of an open quantum system.

# Arguments

  - `H::QuantumObject`: Hamiltonian of the system.
  - `ψ0::QuantumObject`: Initial state of the system.
  - `t_l::AbstractVector`: List of times at which to save the state of the system.
  - `c_ops::Vector`: List of collapse operators.
  - `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
  - `e_ops::Vector`: List of operators for which to calculate expectation values.
  - `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
  - `params::NamedTuple`: Dictionary of parameters to pass to the solver.
  - `seeds::Union{Nothing, Vector{Int}}`: List of seeds for the random number generator. Length must be equal to the number of trajectories provided.
  - `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
  - `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns

  - `prob::ODEProblem`: The ODEProblem for the Monte Carlo wave function time evolution.
"""
function mcsolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{Te,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    seeds::Union{Nothing,Vector{Int}} = nothing,
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix,Te<:AbstractMatrix,TJC<:LindbladJumpCallbackType}
    H_eff = H - T2(0.5im) * mapreduce(op -> op' * op, +, c_ops)

    e_ops2 = Vector{Te}(undef, length(e_ops))
    for i in eachindex(e_ops)
        e_ops2[i] = get_data(e_ops[i])
    end

    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    cache_mc = similar(ψ0.data)
    weights_mc = Array{Float64}(undef, length(c_ops))
    cumsum_weights_mc = similar(weights_mc)

    jump_times_which_init_size = 200
    jump_times = Vector{Float64}(undef, jump_times_which_init_size)
    jump_which = Vector{Int16}(undef, jump_times_which_init_size)

    params2 = (
        expvals = expvals,
        e_ops_mc = e_ops2,
        is_empty_e_ops_mc = isempty(e_ops),
        progr_mc = ProgressBar(length(t_l), enable = false),
        seeds = seeds,
        random_n = Ref(rand()),
        c_ops = get_data.(c_ops),
        cache_mc = cache_mc,
        weights_mc = weights_mc,
        cumsum_weights_mc = cumsum_weights_mc,
        jump_times = jump_times,
        jump_which = jump_which,
        jump_times_which_init_size = jump_times_which_init_size,
        jump_times_which_idx = Ref(1),
        params...,
    )

    return mcsolveProblem(H_eff, ψ0, t_l, alg, H_t, params2, jump_callback; kwargs...)
end

function mcsolveProblem(
    H_eff::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum},
    params::NamedTuple,
    jump_callback::DiscreteLindbladJumpCallback;
    kwargs...,
) where {T1,T2}
    cb1 = DiscreteCallback(LindbladJumpDiscreteCondition, LindbladJumpAffect!, save_positions = (false, false))
    cb2 = PresetTimeCallback(t_l, _save_func_mcsolve, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, cb2, kwargs2.callback),)) :
        merge(kwargs2, (callback = CallbackSet(cb1, cb2),))

    return sesolveProblem(H_eff, ψ0, t_l; alg = alg, H_t = H_t, params = params, kwargs2...)
end

function mcsolveProblem(
    H_eff::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum},
    params::NamedTuple,
    jump_callback::ContinuousLindbladJumpCallback;
    kwargs...,
) where {T1,T2}
    cb1 = ContinuousCallback(
        LindbladJumpContinuousCondition,
        LindbladJumpAffect!,
        nothing,
        interp_points = jump_callback.interp_points,
        save_positions = (false, false),
    )
    cb2 = PresetTimeCallback(t_l, _save_func_mcsolve, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, cb2, kwargs2.callback),)) :
        merge(kwargs2, (callback = CallbackSet(cb1, cb2),))

    return sesolveProblem(H_eff, ψ0, t_l; alg = alg, H_t = H_t, params = params, kwargs2...)
end

"""
    mcsolveEnsembleProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        t_l::AbstractVector,
        c_ops::Vector{QuantumObject{Tc, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[];
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Vector{QuantumObject{Te, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        jump_callback::TJC=ContinuousLindbladJumpCallback(),
        prob_func::Function=_mcsolve_prob_func,
        output_func::Function=_mcsolve_output_func,
        kwargs...)

Generates the ODEProblem for an ensemble of trajectories of the Monte Carlo wave function
time evolution of an open quantum system.

# Arguments

  - `H::QuantumObject`: Hamiltonian of the system.
  - `ψ0::QuantumObject`: Initial state of the system.
  - `t_l::AbstractVector`: List of times at which to save the state of the system.
  - `c_ops::Vector`: List of collapse operators.
  - `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
  - `e_ops::Vector`: List of operators for which to calculate expectation values.
  - `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
  - `params::NamedTuple`: Dictionary of parameters to pass to the solver.
  - `seeds::Union{Nothing, Vector{Int}}`: List of seeds for the random number generator. Length must be equal to the number of trajectories provided.
  - `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
  - `prob_func::Function`: Function to use for generating the ODEProblem.
  - `output_func::Function`: Function to use for generating the output of a single trajectory.
  - `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns

  - `prob::EnsembleProblem with ODEProblem`: The Ensemble ODEProblem for the Monte Carlo
    wave function time evolution.
"""
function mcsolveEnsembleProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{Te,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    seeds::Union{Nothing,Vector{Int}} = nothing,
    prob_func::Function = _mcsolve_prob_func,
    output_func::Function = _mcsolve_output_func,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix,Te<:AbstractMatrix,TJC<:LindbladJumpCallbackType}
    prob_mc = mcsolveProblem(
        H,
        ψ0,
        t_l,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        seeds = seeds,
        jump_callback = jump_callback,
        kwargs...,
    )

    ensemble_prob = EnsembleProblem(prob_mc, prob_func = prob_func, output_func = output_func, safetycopy = false)

    return ensemble_prob
end

"""
    mcsolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        t_l::AbstractVector,
        c_ops::Vector{QuantumObject{Tc, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[];
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Vector{QuantumObject{Te, OperatorQuantumObject}}=QuantumObject{Matrix, OperatorQuantumObject}[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        n_traj::Int=1,
        ensemble_method=EnsembleThreads(),
        jump_callback::TJC=ContinuousLindbladJumpCallback(),
        kwargs...)

Time evolution of an open quantum system using quantum trajectories.

# Arguments

  - `H::QuantumObject`: Hamiltonian of the system.
  - `ψ0::QuantumObject`: Initial state of the system.
  - `t_l::AbstractVector`: List of times at which to save the state of the system.
  - `c_ops::Vector`: List of collapse operators.
  - `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
  - `e_ops::Vector`: List of operators for which to calculate expectation values.
  - `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
  - `params::NamedTuple`: Dictionary of parameters to pass to the solver.
  - `seeds::Union{Nothing, Vector{Int}}`: List of seeds for the random number generator. Length must be equal to the number of trajectories provided.
  - `n_traj::Int`: Number of trajectories to use.
  - `ensemble_method`: Ensemble method to use.
  - `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
  - `prob_func::Function`: Function to use for generating the ODEProblem.
  - `output_func::Function`: Function to use for generating the output of a single trajectory.
  - `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns

  - `sol::TimeEvolutionMCSol`: The solution of the time evolution.

# Notes

`ensemble_method` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`.
"""
function mcsolve(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{Te,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    seeds::Union{Nothing,Vector{Int}} = nothing,
    n_traj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    prob_func::Function = _mcsolve_prob_func,
    output_func::Function = _mcsolve_output_func,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix,Te<:AbstractMatrix,TJC<:LindbladJumpCallbackType}
    if !isnothing(seeds) && length(seeds) != n_traj
        throw(ArgumentError("Length of seeds must match n_traj ($n_traj), but got $(length(seeds))"))
    end

    ens_prob_mc = mcsolveEnsembleProblem(
        H,
        ψ0,
        t_l,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        seeds = seeds,
        jump_callback = jump_callback,
        prob_func = prob_func,
        output_func = output_func,
        kwargs...,
    )

    return mcsolve(ens_prob_mc; alg = alg, n_traj = n_traj, ensemble_method = ensemble_method, kwargs...)
end

function mcsolve(
    ens_prob_mc::EnsembleProblem;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    n_traj::Int = 1,
    ensemble_method = EnsembleThreads(),
    kwargs...,
)
    sol = solve(ens_prob_mc, alg, ensemble_method, trajectories = n_traj)

    expvals_all = Array{ComplexF64}(undef, length(sol), size(sol[:, 1].prob.p.expvals)...)
    times = Vector{Vector{Float64}}(undef, length(sol))
    states =
        haskey(sol[:, 1].prob.kwargs, :save_idxs) ? Vector{Vector{eltype(sol[:, 1].u[1])}}(undef, length(sol)) :
        Vector{Vector{QuantumObject}}(undef, length(sol))
    jump_times = Vector{Vector{Float64}}(undef, length(sol))
    jump_which = Vector{Vector{Int16}}(undef, length(sol))

    foreach(
        i -> _mcsolve_generate_statistics(sol, i, times, states, expvals_all, jump_times, jump_which),
        eachindex(sol),
    )
    expvals = dropdims(sum(expvals_all, dims = 1), dims = 1) ./ length(sol)

    return TimeEvolutionMCSol(times, states, expvals, expvals_all, jump_times, jump_which)
end
