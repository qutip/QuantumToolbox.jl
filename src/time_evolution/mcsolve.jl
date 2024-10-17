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
    c_ops_herm = internal_params.c_ops_herm
    cache_mc = internal_params.cache_mc
    weights_mc = internal_params.weights_mc
    cumsum_weights_mc = internal_params.cumsum_weights_mc
    random_n = internal_params.random_n
    jump_times = internal_params.jump_times
    jump_which = internal_params.jump_which
    traj_rng = internal_params.traj_rng
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        weights_mc[i] = real(dot(ψ, c_ops_herm[i], ψ))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    r = rand(traj_rng) * sum(weights_mc)
    collaps_idx = getindex(1:length(weights_mc), findfirst(>(r), cumsum_weights_mc))
    mul!(cache_mc, c_ops[collaps_idx], ψ)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)

    random_n[] = rand(traj_rng)
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

    global_rng = internal_params.global_rng
    seed = internal_params.seeds[i]
    traj_rng = typeof(global_rng)()
    seed!(traj_rng, seed)

    prm = merge(
        internal_params,
        (
            expvals = similar(internal_params.expvals),
            cache_mc = similar(internal_params.cache_mc),
            weights_mc = similar(internal_params.weights_mc),
            cumsum_weights_mc = similar(internal_params.weights_mc),
            traj_rng = traj_rng,
            random_n = Ref(rand(traj_rng)),
            progr_mc = ProgressBar(size(internal_params.expvals, 2), enable = false),
            jump_times_which_idx = Ref(1),
            jump_times = similar(internal_params.jump_times),
            jump_which = similar(internal_params.jump_which),
        ),
    )

    return remake(prob, p = prm)
end

# Standard output function
function _mcsolve_output_func(sol, i)
    resize!(sol.prob.p.jump_times, sol.prob.p.jump_times_which_idx[] - 1)
    resize!(sol.prob.p.jump_which, sol.prob.p.jump_times_which_idx[] - 1)
    return (sol, false)
end

# Output function with progress bar update
function _mcsolve_output_func_progress(sol, i)
    next!(sol.prob.p.progr_trajectories)
    return _mcsolve_output_func(sol, i)
end

# Output function with distributed channel update for progress bar
function _mcsolve_output_func_distributed(sol, i)
    put!(sol.prob.p.progr_channel, true)
    return _mcsolve_output_func(sol, i)
end

_mcsolve_dispatch_output_func() = _mcsolve_output_func
_mcsolve_dispatch_output_func(::ET) where {ET<:Union{EnsembleSerial,EnsembleThreads}} = _mcsolve_output_func_progress
_mcsolve_dispatch_output_func(::EnsembleDistributed) = _mcsolve_output_func_distributed

function _mcsolve_generate_statistics(sol, i, states, expvals_all, jump_times, jump_which)
    sol_i = sol[:, i]
    !isempty(sol_i.prob.kwargs[:saveat]) ?
    states[i] = [QuantumObject(normalize!(sol_i.u[i]), dims = sol_i.prob.p.Hdims) for i in 1:length(sol_i.u)] : nothing

    copyto!(view(expvals_all, i, :, :), sol_i.prob.p.expvals)
    jump_times[i] = sol_i.prob.p.jump_times
    return jump_which[i] = sol_i.prob.p.jump_which
end

_mcsolve_make_Heff_QobjEvo(H::QuantumObject, c_ops) = QobjEvo(H - 1im * mapreduce(op -> op' * op, +, c_ops) / 2)
_mcsolve_make_Heff_QobjEvo(H::Tuple, c_ops) = QobjEvo((H..., -1im * mapreduce(op -> op' * op, +, c_ops) / 2))
_mcsolve_make_Heff_QobjEvo(H::QuantumObjectEvolution, c_ops) =
    H + QobjEvo(mapreduce(op -> op' * op, +, c_ops), -1im / 2)

@doc raw"""
    mcsolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple}=nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        rng::AbstractRNG=default_rng(),
        jump_callback::TJC=ContinuousLindbladJumpCallback(),
        kwargs...)

Generates the ODEProblem for a single trajectory of the Monte Carlo wave function time evolution of an open quantum system.

Given a system Hamiltonian ``\hat{H}`` and a list of collapse (jump) operators ``\{\hat{C}_n\}_n``, the evolution of the state ``|\psi(t)\rangle`` is governed by the Schrodinger equation:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle= -i \hat{H}_{\textrm{eff}} |\psi(t)\rangle
```

with a non-Hermitian effective Hamiltonian:

```math
\hat{H}_{\textrm{eff}} = \hat{H} - \frac{i}{2} \sum_n \hat{C}_n^\dagger \hat{C}_n.
```

To the first-order of the wave function in a small time ``\delta t``, the strictly negative non-Hermitian portion in ``\hat{H}_{\textrm{eff}}`` gives rise to a reduction in the norm of the wave function, namely

```math
\langle \psi(t+\delta t) | \psi(t+\delta t) \rangle = 1 - \delta p,
```

where 

```math
\delta p = \delta t \sum_n \langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle
```

is the corresponding quantum jump probability.

If the environmental measurements register a quantum jump, the wave function undergoes a jump into a state defined by projecting ``|\psi(t)\rangle`` using the collapse operator ``\hat{C}_n`` corresponding to the measurement, namely

```math
| \psi(t+\delta t) \rangle = \frac{\hat{C}_n |\psi(t)\rangle}{ \sqrt{\langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle} }
```

# Arguments

- `H::QuantumObject`: Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::Union{Nothing,AbstractVector,Tuple}`: List of collapse operators ``\{\hat{C}_n\}_n``.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `rng::AbstractRNG`: Random number generator for reproducibility.
- `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::ODEProblem`: The ODEProblem for the Monte Carlo wave function time evolution.
"""
function mcsolveProblem(
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    kwargs...,
) where {DT1,DT2,TJC<:LindbladJumpCallbackType}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    c_ops isa Nothing &&
        throw(ArgumentError("The list of collapse operators must be provided. Use sesolveProblem instead."))

    tlist = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    H_eff_evo = _mcsolve_make_Heff_QobjEvo(H, c_ops)

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(tlist))
        is_empty_e_ops_mc = true
        e_ops_data = ()
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(tlist))
        e_ops_data = get_data.(e_ops)
        is_empty_e_ops_mc = isempty(e_ops)
    end

    saveat = is_empty_e_ops_mc ? tlist : [tlist[end]]
    # We disable the progress bar of the sesolveProblem because we use a global progress bar for all the trajectories
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat, progress_bar = Val(false))
    kwargs2 = merge(default_values, kwargs)

    cache_mc = similar(ψ0.data)
    weights_mc = Array{Float64}(undef, length(c_ops))
    cumsum_weights_mc = similar(weights_mc)

    jump_times_which_init_size = 200
    jump_times = Vector{Float64}(undef, jump_times_which_init_size)
    jump_which = Vector{Int16}(undef, jump_times_which_init_size)

    c_ops_data = get_data.(c_ops)
    c_ops_herm_data = map(op -> op' * op, c_ops_data)

    params2 = (
        expvals = expvals,
        e_ops_mc = e_ops_data,
        is_empty_e_ops_mc = is_empty_e_ops_mc,
        progr_mc = ProgressBar(length(tlist), enable = false),
        traj_rng = rng,
        c_ops = c_ops_data,
        c_ops_herm = c_ops_herm_data,
        cache_mc = cache_mc,
        weights_mc = weights_mc,
        cumsum_weights_mc = cumsum_weights_mc,
        jump_times = jump_times,
        jump_which = jump_which,
        jump_times_which_init_size = jump_times_which_init_size,
        jump_times_which_idx = Ref(1),
        params...,
    )

    return mcsolveProblem(H_eff_evo, ψ0, tlist, alg, params2, jump_callback; kwargs2...)
end

function mcsolveProblem(
    H_eff_evo::QuantumObjectEvolution{DT1,OperatorQuantumObject},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    alg::OrdinaryDiffEqAlgorithm,
    params::NamedTuple,
    jump_callback::DiscreteLindbladJumpCallback;
    kwargs...,
) where {DT1,DT2}
    cb1 = DiscreteCallback(LindbladJumpDiscreteCondition, LindbladJumpAffect!, save_positions = (false, false))
    cb2 = PresetTimeCallback(tlist, _save_func_mcsolve, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, cb2, kwargs2.callback),)) :
        merge(kwargs2, (callback = CallbackSet(cb1, cb2),))

    return sesolveProblem(H_eff_evo, ψ0, tlist; alg = alg, params = params, kwargs2...)
end

function mcsolveProblem(
    H_eff_evo::QuantumObjectEvolution{DT1,OperatorQuantumObject},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    alg::OrdinaryDiffEqAlgorithm,
    params::NamedTuple,
    jump_callback::ContinuousLindbladJumpCallback;
    kwargs...,
) where {DT1,DT2}
    cb1 = ContinuousCallback(
        LindbladJumpContinuousCondition,
        LindbladJumpAffect!,
        nothing,
        interp_points = jump_callback.interp_points,
        save_positions = (false, false),
    )
    cb2 = PresetTimeCallback(tlist, _save_func_mcsolve, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb1, cb2, kwargs2.callback),)) :
        merge(kwargs2, (callback = CallbackSet(cb1, cb2),))

    return sesolveProblem(H_eff_evo, ψ0, tlist; alg = alg, params = params, kwargs2...)
end

@doc raw"""
    mcsolveEnsembleProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple}=nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        ntraj::Int=1,
        ensemble_method=EnsembleThreads(),
        jump_callback::TJC=ContinuousLindbladJumpCallback(),
        prob_func::Function=_mcsolve_prob_func,
        output_func::Function=_mcsolve_output_func,
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Generates the `EnsembleProblem` of `ODEProblem`s for the ensemble of trajectories of the Monte Carlo wave function time evolution of an open quantum system.

Given a system Hamiltonian ``\hat{H}`` and a list of collapse (jump) operators ``\{\hat{C}_n\}_n``, the evolution of the state ``|\psi(t)\rangle`` is governed by the Schrodinger equation:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle= -i \hat{H}_{\textrm{eff}} |\psi(t)\rangle
```

with a non-Hermitian effective Hamiltonian:

```math
\hat{H}_{\textrm{eff}} = \hat{H} - \frac{i}{2} \sum_n \hat{C}_n^\dagger \hat{C}_n.
```

To the first-order of the wave function in a small time ``\delta t``, the strictly negative non-Hermitian portion in ``\hat{H}_{\textrm{eff}}`` gives rise to a reduction in the norm of the wave function, namely

```math
\langle \psi(t+\delta t) | \psi(t+\delta t) \rangle = 1 - \delta p,
```

where 

```math
\delta p = \delta t \sum_n \langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle
```

is the corresponding quantum jump probability.

If the environmental measurements register a quantum jump, the wave function undergoes a jump into a state defined by projecting ``|\psi(t)\rangle`` using the collapse operator ``\hat{C}_n`` corresponding to the measurement, namely

```math
| \psi(t+\delta t) \rangle = \frac{\hat{C}_n |\psi(t)\rangle}{ \sqrt{\langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle} }
```

# Arguments

- `H::QuantumObject`: Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::Union{Nothing,AbstractVector,Tuple}`: List of collapse operators ``\{\hat{C}_n\}_n``.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `rng::AbstractRNG`: Random number generator for reproducibility.
- `ntraj::Int`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use.
- `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
- `prob_func::Function`: Function to use for generating the ODEProblem.
- `output_func::Function`: Function to use for generating the output of a single trajectory.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::EnsembleProblem with ODEProblem`: The Ensemble ODEProblem for the Monte Carlo wave function time evolution.
"""
function mcsolveEnsembleProblem(
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    prob_func::Function = _mcsolve_prob_func,
    output_func::Function = _mcsolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2,TJC<:LindbladJumpCallbackType}
    progr = ProgressBar(ntraj, enable = getVal(progress_bar))
    if ensemble_method isa EnsembleDistributed
        progr_channel::RemoteChannel{Channel{Bool}} = RemoteChannel(() -> Channel{Bool}(1))
        @async while take!(progr_channel)
            next!(progr)
        end
        params = merge(params, (progr_channel = progr_channel,))
    else
        params = merge(params, (progr_trajectories = progr,))
    end

    # Stop the async task if an error occurs
    try
        seeds = map(i -> rand(rng, UInt64), 1:ntraj)
        prob_mc = mcsolveProblem(
            H,
            ψ0,
            tlist,
            c_ops;
            alg = alg,
            e_ops = e_ops,
            params = merge(params, (global_rng = rng, seeds = seeds)),
            rng = rng,
            jump_callback = jump_callback,
            kwargs...,
        )

        ensemble_prob = EnsembleProblem(prob_mc, prob_func = prob_func, output_func = output_func, safetycopy = false)

        return ensemble_prob
    catch e
        if ensemble_method isa EnsembleDistributed
            put!(progr_channel, false)
        end
        rethrow()
    end
end

@doc raw"""
    mcsolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        prob_func::Function = _mcsolve_prob_func,
        output_func::Function = _mcsolve_dispatch_output_func(ensemble_method),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Time evolution of an open quantum system using quantum trajectories.

Given a system Hamiltonian ``\hat{H}`` and a list of collapse (jump) operators ``\{\hat{C}_n\}_n``, the evolution of the state ``|\psi(t)\rangle`` is governed by the Schrodinger equation:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle= -i \hat{H}_{\textrm{eff}} |\psi(t)\rangle
```

with a non-Hermitian effective Hamiltonian:

```math
\hat{H}_{\textrm{eff}} = \hat{H} - \frac{i}{2} \sum_n \hat{C}_n^\dagger \hat{C}_n.
```

To the first-order of the wave function in a small time ``\delta t``, the strictly negative non-Hermitian portion in ``\hat{H}_{\textrm{eff}}`` gives rise to a reduction in the norm of the wave function, namely

```math
\langle \psi(t+\delta t) | \psi(t+\delta t) \rangle = 1 - \delta p,
```

where 

```math
\delta p = \delta t \sum_n \langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle
```

is the corresponding quantum jump probability.

If the environmental measurements register a quantum jump, the wave function undergoes a jump into a state defined by projecting ``|\psi(t)\rangle`` using the collapse operator ``\hat{C}_n`` corresponding to the measurement, namely

```math
| \psi(t+\delta t) \rangle = \frac{\hat{C}_n |\psi(t)\rangle}{ \sqrt{\langle \psi(t) | \hat{C}_n^\dagger \hat{C}_n | \psi(t) \rangle} }
```

# Arguments

- `H::QuantumObject`: Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::Union{Nothing,AbstractVector,Tuple}`: List of collapse operators ``\{\hat{C}_n\}_n``.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `rng::AbstractRNG`: Random number generator for reproducibility.
- `ntraj::Int`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use.
- `jump_callback::LindbladJumpCallbackType`: The Jump Callback type: Discrete or Continuous.
- `prob_func::Function`: Function to use for generating the ODEProblem.
- `output_func::Function`: Function to use for generating the output of a single trajectory.
- `kwargs...`: Additional keyword arguments to pass to the solver.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.

# Notes

- `ensemble_method` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`
- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionMCSol`: The solution of the time evolution. See also [`TimeEvolutionMCSol`](@ref)
"""
function mcsolve(
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    prob_func::Function = _mcsolve_prob_func,
    output_func::Function = _mcsolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2,TJC<:LindbladJumpCallbackType}
    ens_prob_mc = mcsolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        params = params,
        rng = rng,
        ntraj = ntraj,
        ensemble_method = ensemble_method,
        jump_callback = jump_callback,
        prob_func = prob_func,
        output_func = output_func,
        progress_bar = progress_bar,
        kwargs...,
    )

    return mcsolve(ens_prob_mc, tlist; alg = alg, ntraj = ntraj, ensemble_method = ensemble_method)
end

function mcsolve(
    ens_prob_mc::EnsembleProblem,
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
)
    try
        sol = solve(ens_prob_mc, alg, ensemble_method, trajectories = ntraj)

        if ensemble_method isa EnsembleDistributed
            put!(sol[:, 1].prob.p.progr_channel, false)
        end

        _sol_1 = sol[:, 1]

        expvals_all = Array{ComplexF64}(undef, length(sol), size(_sol_1.prob.p.expvals)...)
        states =
            isempty(_sol_1.prob.kwargs[:saveat]) ? fill(QuantumObject[], length(sol)) :
            Vector{Vector{QuantumObject}}(undef, length(sol))
        jump_times = Vector{Vector{Float64}}(undef, length(sol))
        jump_which = Vector{Vector{Int16}}(undef, length(sol))

        foreach(i -> _mcsolve_generate_statistics(sol, i, states, expvals_all, jump_times, jump_which), eachindex(sol))
        expvals = dropdims(sum(expvals_all, dims = 1), dims = 1) ./ length(sol)

        return TimeEvolutionMCSol(
            ntraj,
            tlist,
            states,
            expvals,
            expvals_all,
            jump_times,
            jump_which,
            sol.converged,
            _sol_1.alg,
            _sol_1.prob.kwargs[:abstol],
            _sol_1.prob.kwargs[:reltol],
        )
    catch e
        if ensemble_method isa EnsembleDistributed
            put!(ens_prob_mc.prob.p.progr_channel, false)
        end
        rethrow()
    end
end
