export mcsolveProblem, mcsolveEnsembleProblem, mcsolve
export ContinuousLindbladJumpCallback, DiscreteLindbladJumpCallback

function _mcsolve_prob_func(prob, i, repeat, global_rng, seeds, tlist)
    seed = seeds[i]
    traj_rng = typeof(global_rng)()
    seed!(traj_rng, seed)

    f = deepcopy(prob.f.f)
    cb = _mcsolve_initialize_callbacks(prob, tlist, traj_rng)

    return remake(prob, f = f, callback = cb)
end

# Standard output function
function _mcsolve_output_func(sol, i)
    idx = _mc_get_jump_callback(sol).affect!.col_times_which_idx[]
    resize!(_mc_get_jump_callback(sol).affect!.col_times, idx - 1)
    resize!(_mc_get_jump_callback(sol).affect!.col_which, idx - 1)
    return (sol, false)
end

function _normalize_state!(u, dims, normalize_states)
    getVal(normalize_states) && normalize!(u)
    return QuantumObject(u, Ket(), dims)
end

function _mcsolve_make_Heff_QobjEvo(H::QuantumObject, c_ops)
    c_ops isa Nothing && return QobjEvo(H)
    return QobjEvo(H - 1im * mapreduce(op -> op' * op, +, c_ops) / 2)
end
function _mcsolve_make_Heff_QobjEvo(H::Tuple, c_ops)
    c_ops isa Nothing && return QobjEvo(H)
    return QobjEvo((H..., -1im * mapreduce(op -> op' * op, +, c_ops) / 2))
end
function _mcsolve_make_Heff_QobjEvo(H::QuantumObjectEvolution, c_ops)
    c_ops isa Nothing && return H
    return H + QobjEvo(mapreduce(op -> op' * op, +, c_ops), -1im / 2)
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
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
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
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    kwargs...,
) where {TJC<:LindbladJumpCallbackType}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    c_ops isa Nothing &&
        throw(ArgumentError("The list of collapse operators must be provided. Use sesolveProblem instead."))

    tlist = _check_tlist(tlist, _float_type(ψ0))

    H_eff_evo = _mcsolve_make_Heff_QobjEvo(H, c_ops)

    T = Base.promote_eltype(H_eff_evo, ψ0)

    # We disable the progress bar of the sesolveProblem because we use a global progress bar for all the trajectories
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., progress_bar = Val(false))
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
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemblealg`: Ensemble algorithm to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
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
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    kwargs...,
) where {TJC<:LindbladJumpCallbackType}
    _prob_func = isnothing(prob_func) ? _ensemble_dispatch_prob_func(rng, ntraj, tlist, _mcsolve_prob_func) : prob_func
    _output_func =
        output_func isa Nothing ?
        _ensemble_dispatch_output_func(ensemblealg, progress_bar, ntraj, _mcsolve_output_func) : output_func

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
        prob_mc.dimensions,
        (progr = _output_func[2], channel = _output_func[3]),
    )

    return ensemble_prob
end

@doc raw"""
    mcsolve(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
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
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm to use for the ODE solver. Default to `Tsit5()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemblealg`: Ensemble algorithm to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
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
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 500,
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    progress_bar::Union{Val,Bool} = Val(true),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    keep_runs_results::Union{Val,Bool} = Val(false),
    normalize_states::Union{Val,Bool} = Val(true),
    kwargs...,
) where {TJC<:LindbladJumpCallbackType}
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
        ensemblealg = ensemblealg,
        jump_callback = jump_callback,
        progress_bar = progress_bar,
        prob_func = prob_func,
        output_func = output_func,
        kwargs...,
    )

    return mcsolve(ens_prob_mc, alg, ntraj, ensemblealg, makeVal(keep_runs_results), normalize_states)
end

function mcsolve(
    ens_prob_mc::TimeEvolutionProblem,
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    ntraj::Int = 500,
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    keep_runs_results = Val(false),
    normalize_states = Val(true),
)
    sol = _ensemble_dispatch_solve(ens_prob_mc, alg, ensemblealg, ntraj)

    dims = ens_prob_mc.dimensions
    _sol_1 = sol[:, 1]
    _expvals_sol_1 = _get_expvals(_sol_1, SaveFuncMCSolve)

    _expvals_all =
        _expvals_sol_1 isa Nothing ? nothing : map(i -> _get_expvals(sol[:, i], SaveFuncMCSolve), eachindex(sol))
    expvals_all = _expvals_all isa Nothing ? nothing : stack(_expvals_all, dims = 2) # Stack on dimension 2 to align with QuTiP

    # stack to transform Vector{Vector{QuantumObject}} -> Matrix{QuantumObject}
    states_all = stack(map(i -> _normalize_state!.(sol[:, i].u, Ref(dims), normalize_states), eachindex(sol)), dims = 1)

    col_times = map(i -> _mc_get_jump_callback(sol[:, i]).affect!.col_times, eachindex(sol))
    col_which = map(i -> _mc_get_jump_callback(sol[:, i]).affect!.col_which, eachindex(sol))

    kwargs = NamedTuple(_sol_1.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionMCSol(
        ntraj,
        ens_prob_mc.times,
        _sol_1.t,
        _store_multitraj_states(states_all, keep_runs_results),
        _store_multitraj_expect(expvals_all, keep_runs_results),
        col_times,
        col_which,
        sol.converged,
        _sol_1.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end
