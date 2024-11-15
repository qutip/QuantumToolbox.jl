export mcsolveProblem, mcsolveEnsembleProblem, mcsolve
export ContinuousLindbladJumpCallback, DiscreteLindbladJumpCallback

const jump_times_which_init_size = 200

function _save_func_mcsolve(integrator, e_ops, is_empty_e_ops)
    expvals = integrator.p.expvals
    progr = integrator.p.progr
    cache_mc = integrator.p.mcsolve_params.cache_mc
    if !is_empty_e_ops
        copyto!(cache_mc, integrator.u)
        normalize!(cache_mc)
        ψ = cache_mc
        _expect = op -> dot(ψ, get_data(op), ψ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    next!(progr)
    u_modified!(integrator, false)
    return nothing
end

function LindbladJumpAffect!(integrator, c_ops, c_ops_herm)
    params = integrator.p
    cache_mc = params.mcsolve_params.cache_mc
    weights_mc = params.mcsolve_params.weights_mc
    cumsum_weights_mc = params.mcsolve_params.cumsum_weights_mc
    random_n = params.mcsolve_params.random_n
    jump_times = params.mcsolve_params.jump_times
    jump_which = params.mcsolve_params.jump_which
    jump_times_which_idx = params.mcsolve_params.jump_times_which_idx
    traj_rng = params.mcsolve_params.traj_rng
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        weights_mc[i] = real(dot(ψ, c_ops_herm[i], ψ))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    r = rand(traj_rng) * sum(real, weights_mc)
    collapse_idx = getindex(1:length(weights_mc), findfirst(x -> real(x) > r, cumsum_weights_mc))
    mul!(cache_mc, c_ops[collapse_idx], ψ)
    normalize!(cache_mc)
    copyto!(integrator.u, cache_mc)

    @inbounds random_n[1] = rand(traj_rng)

    @inbounds idx = round(Int, real(jump_times_which_idx[1]))
    @inbounds jump_times[idx] = integrator.t
    @inbounds jump_which[idx] = collapse_idx
    @inbounds jump_times_which_idx[1] += 1
    @inbounds if real(jump_times_which_idx[1]) > length(jump_times)
        resize!(jump_times, length(jump_times) + jump_times_which_init_size)
        resize!(jump_which, length(jump_which) + jump_times_which_init_size)
    end
end

_mcsolve_continuous_condition(u, t, integrator) =
    @inbounds real(integrator.p.mcsolve_params.random_n[1]) - real(dot(u, u))

_mcsolve_discrete_condition(u, t, integrator) =
    @inbounds real(dot(u, u)) < real(integrator.p.mcsolve_params.random_n[1])

function _mcsolve_prob_func(prob, i, repeat, global_rng, seeds)
    params = prob.p

    seed = seeds[i]
    traj_rng = typeof(global_rng)()
    seed!(traj_rng, seed)

    expvals = similar(params.expvals)
    progr = ProgressBar(size(expvals, 2), enable = false)

    T = eltype(expvals)

    mcsolve_params = (
        traj_rng = traj_rng,
        random_n = T[rand(traj_rng)],
        cache_mc = similar(params.mcsolve_params.cache_mc),
        weights_mc = similar(params.mcsolve_params.weights_mc),
        cumsum_weights_mc = similar(params.mcsolve_params.weights_mc),
        jump_times = similar(params.mcsolve_params.jump_times),
        jump_which = similar(params.mcsolve_params.jump_which),
        jump_times_which_idx = T[1],
    )

    p = TimeEvolutionParameters(params.params, expvals, progr, mcsolve_params)

    f = deepcopy(prob.f.f)

    return remake(prob, f = f, p = p)
end

function _mcsolve_dispatch_prob_func(rng, ntraj)
    seeds = map(i -> rand(rng, UInt64), 1:ntraj)
    return (prob, i, repeat) -> _mcsolve_prob_func(prob, i, repeat, rng, seeds)
end

# Standard output function
function _mcsolve_output_func(sol, i)
    @inbounds idx = round(Int, real(sol.prob.p.mcsolve_params.jump_times_which_idx[1]))
    resize!(sol.prob.p.mcsolve_params.jump_times, idx - 1)
    resize!(sol.prob.p.mcsolve_params.jump_which, idx - 1)
    return (sol, false)
end

# Output function with progress bar update
function _mcsolve_output_func_progress(sol, i, progr)
    next!(progr)
    return _mcsolve_output_func(sol, i)
end

# Output function with distributed channel update for progress bar
function _mcsolve_output_func_distributed(sol, i, channel)
    put!(channel, true)
    return _mcsolve_output_func(sol, i)
end

function _mcsolve_dispatch_output_func(::ET, progress_bar, ntraj) where {ET<:Union{EnsembleSerial,EnsembleThreads}}
    if getVal(progress_bar)
        progr = ProgressBar(ntraj, enable = getVal(progress_bar))
        f = (sol, i) -> _mcsolve_output_func_progress(sol, i, progr)
        return (f, progr, nothing)
    else
        return (_mcsolve_output_func, nothing, nothing)
    end
end
function _mcsolve_dispatch_output_func(
    ::ET,
    progress_bar,
    ntraj,
) where {ET<:Union{EnsembleSplitThreads,EnsembleDistributed}}
    if getVal(progress_bar)
        progr = ProgressBar(ntraj, enable = getVal(progress_bar))
        progr_channel::RemoteChannel{Channel{Bool}} = RemoteChannel(() -> Channel{Bool}(1))

        f = (sol, i) -> _mcsolve_output_func_distributed(sol, i, progr_channel)
        return (f, progr, progr_channel)
    else
        return (_mcsolve_output_func, nothing, nothing)
    end
end

function _normalize_state!(u, dims, normalize_states)
    getVal(normalize_states) && normalize!(u)
    return QuantumObject(u, type = Ket, dims = dims)
end

function _generate_mcsolve_kwargs(e_ops, tlist, c_ops, jump_callback, kwargs)
    c_ops_data = get_data.(c_ops)
    c_ops_herm_data = map(op -> op' * op, c_ops_data)

    _affect = integrator -> LindbladJumpAffect!(integrator, c_ops_data, c_ops_herm_data)

    if jump_callback isa DiscreteLindbladJumpCallback
        cb1 = DiscreteCallback(_mcsolve_discrete_condition, _affect, save_positions = (false, false))
    else
        cb1 = ContinuousCallback(
            _mcsolve_continuous_condition,
            _affect,
            nothing,
            interp_points = jump_callback.interp_points,
            save_positions = (false, false),
        )
    end

    if e_ops isa Nothing
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb1, kwargs.callback),)) :
            merge(kwargs, (callback = cb1,))
        return kwargs2
    else
        is_empty_e_ops = isempty(e_ops)
        f = integrator -> _save_func_mcsolve(integrator, e_ops, is_empty_e_ops)
        cb2 = PresetTimeCallback(tlist, f, save_positions = (false, false))
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(cb1, cb2, kwargs.callback),)) :
            merge(kwargs, (callback = CallbackSet(cb1, cb2),))
        return kwargs2
    end
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
        H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{DT2,KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
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
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` or `AbstractVector` of parameters to pass to the solver.
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
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
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

    T = Base.promote_eltype(H_eff_evo, ψ0)

    if e_ops isa Nothing
        expvals = Array{T}(undef, 0, length(tlist))
        is_empty_e_ops = true
    else
        expvals = Array{T}(undef, length(e_ops), length(tlist))
        is_empty_e_ops = isempty(e_ops)
    end

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    # We disable the progress bar of the sesolveProblem because we use a global progress bar for all the trajectories
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat, progress_bar = Val(false))
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_mcsolve_kwargs(e_ops, tlist, c_ops, jump_callback, kwargs2)

    cache_mc = similar(ψ0.data, T)
    weights_mc = similar(ψ0.data, T, length(c_ops)) # It should be a Float64 Vector, but we have to keep the same type for all the parameters due to SciMLStructures.jl
    cumsum_weights_mc = similar(weights_mc)

    jump_times = similar(ψ0.data, T, jump_times_which_init_size)
    jump_which = similar(ψ0.data, T, jump_times_which_init_size)
    jump_times_which_idx = T[1] # We could use a Ref, but we have to keep the same type for all the parameters due to SciMLStructures.jl

    random_n = similar(ψ0.data, T, 1) # We could use a Ref, but we have to keep the same type for all the parameters due to SciMLStructures.jl.
    random_n[1] = rand(rng)

    progr = ProgressBar(length(tlist), enable = false)

    mcsolve_params = (
        traj_rng = rng,
        random_n = random_n,
        cache_mc = cache_mc,
        weights_mc = weights_mc,
        cumsum_weights_mc = cumsum_weights_mc,
        jump_times = jump_times,
        jump_which = jump_which,
        jump_times_which_idx = jump_times_which_idx,
    )
    p = TimeEvolutionParameters(params, expvals, progr, mcsolve_params)

    return sesolveProblem(H_eff_evo, ψ0, tlist; params = p, kwargs3...)
end

@doc raw"""
    mcsolveEnsembleProblem(
        H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{DT2,KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
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
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` or `AbstractVector` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
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
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    progress_bar::Union{Val,Bool} = Val(true),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    kwargs...,
) where {DT1,DT2,TJC<:LindbladJumpCallbackType}
    _prob_func = prob_func isa Nothing ? _mcsolve_dispatch_prob_func(rng, ntraj) : prob_func
    _output_func =
        output_func isa Nothing ? _mcsolve_dispatch_output_func(ensemble_method, progress_bar, ntraj) : output_func

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
        EnsembleProblem(prob_mc.prob, prob_func = _prob_func, output_func = _output_func[1], safetycopy = true),
        prob_mc.times,
        prob_mc.dims,
        (progr = _output_func[2], channel = _output_func[3]),
    )

    return ensemble_prob
end

@doc raw"""
    mcsolve(
        H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{DT2,KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        progress_bar::Union{Val,Bool} = Val(true),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
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
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm to use for the ODE solver. Default to `Tsit5()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` or `AbstractVector` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `normalize_states`: Whether to normalize the states. Default to `Val(true)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- `ensemble_method` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`
- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionMCSol`: The solution of the time evolution. See also [`TimeEvolutionMCSol`](@ref).
"""
function mcsolve(
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    progress_bar::Union{Val,Bool} = Val(true),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    normalize_states::Union{Val,Bool} = Val(true),
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
        progress_bar = progress_bar,
        prob_func = prob_func,
        output_func = output_func,
        kwargs...,
    )

    return mcsolve(ens_prob_mc, alg, ntraj, ensemble_method, normalize_states)
end

function _mcsolve_solve_ens(
    ens_prob_mc::TimeEvolutionProblem,
    alg::OrdinaryDiffEqAlgorithm,
    ensemble_method::ET,
    ntraj::Int,
) where {ET<:Union{EnsembleSplitThreads,EnsembleDistributed}}
    sol = nothing

    @sync begin
        @async while take!(ens_prob_mc.kwargs.channel)
            next!(ens_prob_mc.kwargs.progr)
        end

        @async begin
            sol = solve(ens_prob_mc.prob, alg, ensemble_method, trajectories = ntraj)
            put!(ens_prob_mc.kwargs.channel, false)
        end
    end

    return sol
end

function _mcsolve_solve_ens(
    ens_prob_mc::TimeEvolutionProblem,
    alg::OrdinaryDiffEqAlgorithm,
    ensemble_method,
    ntraj::Int,
)
    sol = solve(ens_prob_mc.prob, alg, ensemble_method, trajectories = ntraj)
    return sol
end

function mcsolve(
    ens_prob_mc::TimeEvolutionProblem,
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    normalize_states = Val(true),
)
    sol = _mcsolve_solve_ens(ens_prob_mc, alg, ensemble_method, ntraj)

    dims = ens_prob_mc.dims
    _sol_1 = sol[:, 1]

    expvals_all = mapreduce(i -> sol[:, i].prob.p.expvals, (x, y) -> cat(x, y, dims = 3), eachindex(sol))
    states = map(i -> _normalize_state!.(sol[:, i].u, Ref(dims), normalize_states), eachindex(sol))
    jump_times = map(i -> real.(sol[:, i].prob.p.mcsolve_params.jump_times), eachindex(sol))
    jump_which = map(i -> round.(Int, sol[:, i].prob.p.mcsolve_params.jump_which), eachindex(sol))

    expvals = dropdims(sum(expvals_all, dims = 3), dims = 3) ./ length(sol)

    return TimeEvolutionMCSol(
        ntraj,
        ens_prob_mc.times,
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
end
