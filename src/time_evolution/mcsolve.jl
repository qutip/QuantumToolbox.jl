export mcsolveProblem, mcsolveEnsembleProblem, mcsolve, mcsolve_map
export ContinuousLindbladJumpCallback, DiscreteLindbladJumpCallback

function _mcsolve_prob_func(prob, ctx, tlist; kwargs...)
    f = _copy_for_trajectory(prob.f.f)
    cb = _mcsolve_initialize_callbacks(prob, tlist, ctx.rng)

    return remake(prob, f = f, callback = cb)
end

# Standard output function
function _mcsolve_output_func(sol, ctx)
    idx = _mc_get_jump_callback(sol).affect!.col_times_which_idx[]
    resize!(_mc_get_jump_callback(sol).affect!.col_times, idx - 1)
    resize!(_mc_get_jump_callback(sol).affect!.col_which, idx - 1)
    return (sol, false)
end

# Keep only the data needed to build a `TimeEvolutionMCSol`, dropping the operators and the integrator cache stored in the `ODESolution`.
# This reduces the memory and the data sent back by the workers of distributed ensembles.
function _mcsolve_trajectory_output(sol::AbstractODESolution)
    jump_affect! = _mc_get_jump_callback(sol).affect!
    n_jumps = jump_affect!.col_times_which_idx[] - 1
    return (
        times_states = sol.t,
        states = sol.u,
        expect = _get_expvals(sol, SaveFuncMCSolve),
        col_times = resize!(jump_affect!.col_times, n_jumps),
        col_which = resize!(jump_affect!.col_which, n_jumps),
    )
end

function _normalize_state!(u, dims, normalize_states)
    getVal(normalize_states) && normalize!(u)
    return QuantumObject(u, Ket(), dims)
end

# Build the `TimeEvolutionMCSol` from the outputs of `_mcsolve_trajectory_output` of all the trajectories
function _gen_mcsolve_solution(
        trajs::AbstractVector,
        times,
        dimensions,
        alg,
        abstol,
        reltol,
        converged::Bool,
        normalize_states,
        keep_runs_results,
    )
    traj_1 = first(trajs)

    expvals_all = traj_1.expect isa Nothing ? nothing : stack(map(traj -> traj.expect, trajs), dims = 2) # Stack on dimension 2 to align with QuTiP

    # stack to transform Vector{Vector{QuantumObject}} -> Matrix{QuantumObject}
    states_all = stack(map(traj -> _normalize_state!.(traj.states, Ref(dimensions), normalize_states), trajs), dims = 1)

    return TimeEvolutionMCSol(
        length(trajs),
        times,
        traj_1.times_states,
        _store_multitraj_states(states_all, makeVal(keep_runs_results)),
        _store_multitraj_expect(expvals_all, makeVal(keep_runs_results)),
        map(traj -> traj.col_times, trajs),
        map(traj -> traj.col_which, trajs),
        converged,
        alg,
        abstol,
        reltol,
    )
end

function _mcsolve_make_Heff_QobjEvo(H::QuantumObject, c_ops)
    c_ops isa Nothing && return QuantumObjectEvolution(H)
    return QuantumObjectEvolution(H - 1im * mapreduce(op -> op' * op, +, c_ops) / 2)
end
function _mcsolve_make_Heff_QobjEvo(H::Tuple, c_ops)
    c_ops isa Nothing && return QuantumObjectEvolution(H)
    return QuantumObjectEvolution((H..., -1im * sum(op -> op' * op, c_ops) / 2))
end
function _mcsolve_make_Heff_QobjEvo(H::QuantumObjectEvolution, c_ops)
    c_ops isa Nothing && return H
    return H + QuantumObjectEvolution(sum(op -> -1im * op' * op / 2, c_ops))
end

@doc raw"""
    mcsolveProblem(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        kwargs...,
    )

Generate the ODEProblem for a single trajectory of the Monte Carlo wave function time evolution of an open quantum system.

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

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`. Each element can be a [`QuantumObject`](@ref) or a [`QuantumObjectEvolution`](@ref) (for time-dependent collapse operators).
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `jump_callback`: The Jump Callback type: [`ContinuousLindbladJumpCallback`](@ref) or [`DiscreteLindbladJumpCallback`](@ref). The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the `ODEProblem` for the Monte Carlo wave function time evolution.
"""
function mcsolveProblem(
        H::Union{AbstractQuantumObject{Operator}, Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    c_ops isa Nothing &&
        throw(ArgumentError("The list of collapse operators must be provided. Use sesolveProblem instead."))

    H_eff_evo = _mcsolve_make_Heff_QobjEvo(H, c_ops)

    T = _complex_float_type(Base.promote_eltype(H_eff_evo, ψ0))

    tlist = _check_tlist(tlist, _float_type(T))

    # We disable the progress bar of the sesolveProblem because we use a global progress bar for all the trajectories
    default_values = (default_ode_solver_options(T)..., progress_bar = Val(false))
    kwargs2 = _merge_saveat(tlist, e_ops, default_values; kwargs...)
    kwargs3 = _generate_mcsolve_kwargs(ψ0, T, e_ops, tlist, c_ops, jump_callback, rng, kwargs2)

    return sesolveProblem(H_eff_evo, ψ0, tlist; params = params, kwargs3...)
end

@doc raw"""
    mcsolveEnsembleProblem(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val,Bool} = Val(true),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        kwargs...,
    )

Generate the `EnsembleProblem` of `ODEProblem`s for the ensemble of trajectories of the Monte Carlo wave function time evolution of an open quantum system.

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

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`. Each element can be a [`QuantumObject`](@ref) or a [`QuantumObjectEvolution`](@ref) (for time-dependent collapse operators).
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemblealg`: Ensemble algorithm to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: [`ContinuousLindbladJumpCallback`](@ref) or [`DiscreteLindbladJumpCallback`](@ref). The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `Progress` object, and the (optional) `RemoteChannel` object.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the Ensemble `ODEProblem` for the Monte Carlo wave function time evolution.
"""
function mcsolveEnsembleProblem(
        H::Union{AbstractQuantumObject{Operator}, Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val, Bool} = Val(true),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple, Nothing} = nothing,
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    _prob_func = isnothing(prob_func) ? _ensemble_dispatch_prob_func(tlist, _mcsolve_prob_func) : prob_func
    _output_func =
        output_func isa Nothing ?
        _ensemble_dispatch_output_func(
            ensemblealg,
            progress_bar,
            ntraj,
            _mcsolve_output_func;
            progr_desc = "[mcsolve] ",
        ) : output_func

    prob_mc = mcsolveProblem(
        H,
        ψ0,
        tlist,
        c_ops;
        e_ops = e_ops,
        params = params,
        rng = rng,
        jump_callback = jump_callback,
        kwargs...,
    )

    ensemble_prob = TimeEvolutionProblem(
        EnsembleProblem(prob_mc.prob, prob_func = _prob_func, output_func = _output_func[1], safetycopy = false),
        prob_mc.times,
        prob_mc.states_type,
        prob_mc.dimensions,
        (progr = _output_func[2], channel = _output_func[3], rng = rng, ntraj = ntraj, ensemblealg = ensemblealg),
    )

    return ensemble_prob
end

@doc raw"""
    mcsolve(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val,Bool} = Val(true),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        keep_runs_results::Union{Val,Bool} = Val(false),
        normalize_states::Union{Val,Bool} = Val(true),
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

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`. Each element can be a [`QuantumObject`](@ref) or a [`QuantumObjectEvolution`](@ref) (for time-dependent collapse operators).
- `alg`: The algorithm to use for the ODE solver. Default to `DP5()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemblealg`: Ensemble algorithm to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: [`ContinuousLindbladJumpCallback`](@ref) or [`DiscreteLindbladJumpCallback`](@ref). The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `Progress` object, and the (optional) `RemoteChannel` object.
- `keep_runs_results`: Whether to save the results of each trajectory. Default to `Val(false)`.
- `normalize_states`: Whether to normalize the states. Default to `Val(true)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- `ensemblealg` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`
- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionMCSol`: The solution of the time evolution. See also [`TimeEvolutionMCSol`](@ref).
"""
function mcsolve(
        H::Union{AbstractQuantumObject{Operator}, Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val, Bool} = Val(true),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple, Nothing} = nothing,
        keep_runs_results::Union{Val, Bool} = Val(false),
        normalize_states::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    ens_prob_mc = mcsolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        c_ops;
        e_ops = e_ops,
        params = params,
        rng = rng,
        ntraj = ntraj,
        ensemblealg = ensemblealg,
        jump_callback = jump_callback,
        progress_bar = progress_bar,
        prob_func = prob_func,
        output_func = output_func,
        kwargs...,
    )

    return mcsolve(ens_prob_mc, alg; keep_runs_results = keep_runs_results, normalize_states = normalize_states)
end

function mcsolve(
        ens_prob_mc::TimeEvolutionProblem,
        alg::AbstractODEAlgorithm = DP5();
        keep_runs_results::Union{Val, Bool} = Val(false),
        normalize_states::Union{Val, Bool} = Val(true),
    )
    ntraj = ens_prob_mc.kwargs.ntraj
    sol = _ensemble_dispatch_solve(ens_prob_mc, alg, ens_prob_mc.kwargs.ensemblealg, ntraj; rng = ens_prob_mc.kwargs.rng)

    _sol_1 = sol.u[1]
    kwargs = NamedTuple(_sol_1.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return _gen_mcsolve_solution(
        map(_mcsolve_trajectory_output, sol.u),
        ens_prob_mc.times,
        ens_prob_mc.dimensions,
        _sol_1.alg,
        kwargs.abstol,
        kwargs.reltol,
        sol.converged,
        normalize_states,
        keep_runs_results,
    )
end

@doc raw"""
    mcsolve_map(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::Union{QuantumObject{Ket},AbstractVector{<:QuantumObject{Ket}}},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::AbstractODEAlgorithm = DP5(),
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NullParameters,Tuple} = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val,Bool} = Val(true),
        keep_runs_results::Union{Val,Bool} = Val(false),
        normalize_states::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Solve the quantum trajectories for multiple initial states and parameter sets using ensemble simulation.

This function computes the Monte Carlo wave function time evolution (see [`mcsolve`](@ref)) with `ntraj` trajectories for all combinations (Cartesian product) of initial states and parameter sets. Each trajectory evolves under the non-Hermitian effective Hamiltonian

```math
\hat{H}_{\textrm{eff}} = \hat{H} - \frac{i}{2} \sum_n \hat{C}_n^\dagger \hat{C}_n,
```

interrupted by quantum jumps.

All the trajectories, of all the combinations, are solved within a single `EnsembleProblem`, so that the parallelization (and the load balancing) is performed over the whole set of trajectories at once.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state(s) of the system. Can be a single [`Ket`](@ref) or a `Vector` of [`Ket`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`. Each element can be a [`QuantumObject`](@ref) or a [`QuantumObjectEvolution`](@ref) (for time-dependent or parameter-dependent collapse operators).
- `alg`: The algorithm for the ODE solver. The default is `DP5()`.
- `ensemblealg`: Ensemble algorithm to use for parallel computation. Default is `EnsembleThreads()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: A `Tuple` of parameter sets. Each element should be an `AbstractVector` representing the sweep range for that parameter. The function will solve for all combinations of initial states and parameter sets.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories for each combination of initial state and parameters.
- `jump_callback`: The Jump Callback type: [`ContinuousLindbladJumpCallback`](@ref) or [`DiscreteLindbladJumpCallback`](@ref). The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `keep_runs_results`: Whether to save the results of each trajectory. Default to `Val(false)`.
- `normalize_states`: Whether to normalize the states. Default to `Val(true)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The function returns an array of solutions with dimensions matching the Cartesian product of initial states and parameter sets.
- If `ψ0` is a vector of `m` states and `params = (p1, p2, ...)` where `p1` has length `n1`, `p2` has length `n2`, etc., the output will be of size `(m, n1, n2, ...)`.
- The total number of solved trajectories is `ntraj * m * n1 * n2 * ...`.
- See [`mcsolve`](@ref) for more details.

# Returns

- An array of [`TimeEvolutionMCSol`](@ref) objects with dimensions `(length(ψ0), length(params[1]), length(params[2]), ...)`.
"""
function mcsolve_map(
        H::Union{AbstractQuantumObject{Operator}, Tuple},
        ψ0::AbstractVector{<:QuantumObject{Ket}},
        tlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        alg::AbstractODEAlgorithm = DP5(),
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
        params::Union{NullParameters, Tuple} = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val, Bool} = Val(true),
        keep_runs_results::Union{Val, Bool} = Val(false),
        normalize_states::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    # mapping initial states and parameters
    ψ0_iter = map(state -> to_dense(_complex_float_type(eltype(state)), copy(state.data)), ψ0)
    if params isa NullParameters
        iter = collect(Iterators.product(ψ0_iter, [params])) |> vec # convert nx1 Matrix into Vector
    else
        iter = collect(Iterators.product(ψ0_iter, params...))
    end

    prob = mcsolveProblem(
        H,
        first(ψ0),
        tlist,
        c_ops;
        e_ops = e_ops,
        params = Base.tail(first(iter)),
        rng = rng,
        jump_callback = jump_callback,
        kwargs...,
    )

    return mcsolve_map(
        prob,
        iter,
        alg,
        ensemblealg;
        ntraj = ntraj,
        rng = rng,
        progress_bar = progress_bar,
        keep_runs_results = keep_runs_results,
        normalize_states = normalize_states,
    )
end
mcsolve_map(
    H::Union{AbstractQuantumObject{Operator}, Tuple},
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector,
    c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
    kwargs...,
) = mcsolve_map(H, [ψ0], tlist, c_ops; kwargs...)

# this method is for advanced usage (see `sesolve_map`)
# Each element `(u0, p...)` of `iter` is solved with `ntraj` trajectories.
# A custom `output_func` must return the output of `_mcsolve_trajectory_output`.
function mcsolve_map(
        prob::TimeEvolutionProblem{Ket, <:Dimensions, <:ODEProblem},
        iter::AbstractArray,
        alg::AbstractODEAlgorithm = DP5(),
        ensemblealg::EnsembleAlgorithm = EnsembleThreads();
        ntraj::Int = 500,
        rng::AbstractRNG = default_rng(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple, Nothing} = nothing,
        safetycopy::Union{Bool, Nothing} = nothing,
        progress_bar::Union{Val, Bool} = Val(true),
        keep_runs_results::Union{Val, Bool} = Val(false),
        normalize_states::Union{Val, Bool} = Val(true),
    )
    # The trajectory `sim_id` solves `iter[mod1(sim_id, length(iter))]`, so the order is (A, B, C, A, B, C, ...) instead of (A, A, ..., B, B, ..., C, C, ...).
    # `EnsembleThreads` gives each thread a contiguous block of trajectories, so the expensive elements of `iter` are spread over all the threads.
    ntraj_tot = ntraj * length(iter)
    tlist = prob.times
    _prob_func = isnothing(prob_func) ? (prob, ctx) -> _mcsolve_map_prob_func(prob, ctx, tlist, iter) : prob_func
    _safetycopy = isnothing(safetycopy) ? !isnothing(prob_func) : safetycopy
    _output_func =
        isnothing(output_func) ?
        _ensemble_dispatch_output_func(
            ensemblealg,
            progress_bar,
            ntraj_tot,
            _mcsolve_map_output_func;
            progr_desc = "[mcsolve_map] ",
        ) : output_func
    ens_prob = TimeEvolutionProblem(
        EnsembleProblem(prob.prob, prob_func = _prob_func, output_func = _output_func[1], safetycopy = _safetycopy),
        prob.times,
        prob.states_type,
        prob.dimensions,
        (progr = _output_func[2], channel = _output_func[3]),
    )

    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemblealg, ntraj_tot; rng = rng)

    # handle solution and make it become an Array of TimeEvolutionMCSol
    trajs = reshape(sol.u, length(iter), ntraj) # the i-th row contains the trajectories of iter[i]
    kwargs = NamedTuple(prob.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility
    gen_sol =
        i -> _gen_mcsolve_solution(
        view(trajs, i, :),
        prob.times,
        prob.dimensions,
        alg,
        kwargs.abstol,
        kwargs.reltol,
        sol.converged,
        normalize_states,
        keep_runs_results,
    )

    # We don't use `map` here: `_gen_mcsolve_solution` itself uses `map`, and Julia's inference loses the element type of a `map`
    # whose function calls `map` with a closure (recursion limiting heuristic), making the return type unstable. An explicit loop avoids it.
    sol_1 = gen_sol(1)
    sol_arr = similar(iter, typeof(sol_1))
    sol_arr[1] = sol_1
    for i in 2:length(iter)
        sol_arr[i] = gen_sol(i)
    end
    return sol_arr
end

function _mcsolve_map_prob_func(prob, ctx, tlist, iter)
    x = iter[mod1(ctx.sim_id, length(iter))]
    f = _copy_for_trajectory(prob.f.f)
    cb = _mcsolve_initialize_callbacks(prob, tlist, ctx.rng)

    return remake(prob, f = f, u0 = first(x), p = Base.tail(x), callback = cb)
end

_mcsolve_map_output_func(sol, ctx) = (_mcsolve_trajectory_output(sol), false)
