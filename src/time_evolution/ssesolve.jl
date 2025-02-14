export ssesolveProblem, ssesolveEnsembleProblem, ssesolve

# TODO: Implement the three-argument dot function for SciMLOperators.jl
# Currently, we are assuming a time-independent MatrixOperator
function _ssesolve_update_coeff(u, p, t, op)
    normalize!(u)
    return real(dot(u, op.A, u)) #this is en/2: <Sn + Sn'>/2 = Re<Sn>
end

_ScalarOperator_e(op, f = +) = ScalarOperator(one(eltype(op)), (a, u, p, t) -> f(_ssesolve_update_coeff(u, p, t, op)))

_ScalarOperator_e2_2(op, f = +) =
    ScalarOperator(one(eltype(op)), (a, u, p, t) -> f(_ssesolve_update_coeff(u, p, t, op)^2 / 2))

@doc raw"""
    ssesolveProblem(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        progress_bar::Union{Val,Bool} = Val(true),
        store_measurement::Union{Val, Bool} = Val(false),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{S}_n - \frac{1}{2} \hat{S}_n^\dagger \hat{S}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{S}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{S}_n + \hat{S}_n^\dagger \rangle.
```

Above, ``\hat{S}_n`` are the stochastic collapse operators and ``dW_n(t)`` is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `store_measurement`: Whether to store the measurement results. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The `SDEProblem` for the Stochastic Schrödinger time evolution of the system.
"""
function ssesolveProblem(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    progress_bar::Union{Val,Bool} = Val(true),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
)
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    sc_ops isa Nothing &&
        throw(ArgumentError("The list of stochastic collapse operators must be provided. Use sesolveProblem instead."))

    tlist = _check_tlist(tlist, _FType(ψ0))

    H_eff_evo = _mcsolve_make_Heff_QobjEvo(H, sc_ops)
    isoper(H_eff_evo) || throw(ArgumentError("The Hamiltonian must be an Operator."))
    check_dimensions(H_eff_evo, ψ0)
    dims = H_eff_evo.dimensions

    ψ0 = to_dense(_CType(ψ0), get_data(ψ0))

    progr = ProgressBar(length(tlist), enable = getVal(progress_bar))

    sc_ops_evo_data = Tuple(map(get_data ∘ QobjEvo, sc_ops))

    # Here the coefficients depend on the state, so this is a non-linear operator, which should be implemented with FunctionOperator instead. However, the nonlinearity is only on the coefficients, and it should be safe.
    K_l = sum(
        op -> _ScalarOperator_e(op, +) * op + _ScalarOperator_e2_2(op, -) * IdentityOperator(prod(dims)),
        sc_ops_evo_data,
    )

    K = -1im * get_data(H_eff_evo) + K_l

    D_l = map(op -> op + _ScalarOperator_e(op, -) * IdentityOperator(prod(dims)), sc_ops_evo_data)
    D = DiffusionOperator(D_l)

    kwargs2 = _merge_saveat(tlist, e_ops, DEFAULT_SDE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_stochastic_kwargs(
        e_ops,
        sc_ops,
        makeVal(progress_bar),
        tlist,
        makeVal(store_measurement),
        kwargs2,
        SaveFuncSSESolve,
    )

    tspan = (tlist[1], tlist[end])
    noise = RealWienerProcess!(
        tlist[1],
        zeros(length(sc_ops)),
        zeros(length(sc_ops)),
        save_everystep = getVal(store_measurement),
        rng = rng,
    )
    noise_rate_prototype = similar(ψ0, length(ψ0), length(sc_ops))
    prob = SDEProblem{true}(
        K,
        D,
        ψ0,
        tspan,
        params;
        noise_rate_prototype = noise_rate_prototype,
        noise = noise,
        kwargs3...,
    )

    return TimeEvolutionProblem(prob, tlist, dims)
end

@doc raw"""
    ssesolveEnsembleProblem(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        store_measurement::Union{Val,Bool} = Val(false),
        kwargs...,
    )

Generate the SDE EnsembleProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{S}_n - \frac{1}{2} \hat{S}_n^\dagger \hat{S}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{S}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{S}_n + \hat{S}_n^\dagger \rangle.
```

Above, ``\hat{S}_n`` are the stochastic collapse operators and  ``dW_n(t)`` is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `prob_func`: Function to use for generating the SDEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `store_measurement`: Whether to store the measurement results. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::EnsembleProblem with SDEProblem`: The Ensemble SDEProblem for the Stochastic Shrödinger time evolution.
"""
function ssesolveEnsembleProblem(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
)
    _prob_func =
        isnothing(prob_func) ?
        _ensemble_dispatch_prob_func(
            rng,
            ntraj,
            tlist,
            _stochastic_prob_func;
            n_sc_ops = length(sc_ops),
            store_measurement = makeVal(store_measurement),
        ) : prob_func
    _output_func =
        output_func isa Nothing ?
        _ensemble_dispatch_output_func(ensemble_method, progress_bar, ntraj, _stochastic_output_func) : output_func

    prob_sme = ssesolveProblem(
        H,
        ψ0,
        tlist,
        sc_ops;
        e_ops = e_ops,
        params = params,
        rng = rng,
        progress_bar = Val(false),
        store_measurement = makeVal(store_measurement),
        kwargs...,
    )

    ensemble_prob = TimeEvolutionProblem(
        EnsembleProblem(prob_sme, prob_func = _prob_func, output_func = _output_func[1], safetycopy = true),
        prob_sme.times,
        prob_sme.dimensions,
        (progr = _output_func[2], channel = _output_func[3]),
    )

    return ensemble_prob
end

@doc raw"""
    ssesolve(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::StochasticDiffEqAlgorithm = SRA1(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        store_measurement::Union{Val,Bool} = Val(false),
        kwargs...,
    )


Stochastic Schrödinger equation evolution of a quantum system given the system Hamiltonian ``\hat{H}`` and a list of stochastic collapse (jump) operators ``\{\hat{S}_n\}_n``.
The stochastic evolution of the state ``|\psi(t)\rangle`` is defined by:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{S}_n - \frac{1}{2} \hat{S}_n^\dagger \hat{S}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{S}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{S}_n + \hat{S}_n^\dagger \rangle.
```

Above, ``\hat{S}_n`` are the stochastic collapse operators and ``dW_n(t)`` is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.


# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm to use for the stochastic differential equation. Default is `SRA1()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the SDEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `store_measurement`: Whether to store the measurement results. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (SDE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionStochasticSol`: The solution of the time evolution. See [`TimeEvolutionStochasticSol`](@ref).
"""
function ssesolve(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
)
    ens_prob = ssesolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        sc_ops;
        e_ops = e_ops,
        params = params,
        rng = rng,
        ntraj = ntraj,
        ensemble_method = ensemble_method,
        prob_func = prob_func,
        output_func = output_func,
        progress_bar = progress_bar,
        store_measurement = makeVal(store_measurement),
        kwargs...,
    )

    return ssesolve(ens_prob, alg, ntraj, ensemble_method)
end

function ssesolve(
    ens_prob::TimeEvolutionProblem,
    alg::StochasticDiffEqAlgorithm = SRA1(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
)
    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemble_method, ntraj)

    _sol_1 = sol[:, 1]
    _expvals_sol_1 = _get_expvals(_sol_1, SaveFuncSSESolve)
    _m_expvals_sol_1 = _get_m_expvals(_sol_1, SaveFuncSSESolve)

    normalize_states = Val(false)
    dims = ens_prob.dimensions
    _expvals_all =
        _expvals_sol_1 isa Nothing ? nothing : map(i -> _get_expvals(sol[:, i], SaveFuncSSESolve), eachindex(sol))
    expvals_all = _expvals_all isa Nothing ? nothing : stack(_expvals_all, dims = 2) # Stack on dimension 2 to align with QuTiP
    states = map(i -> _normalize_state!.(sol[:, i].u, Ref(dims), normalize_states), eachindex(sol))

    _m_expvals =
        _m_expvals_sol_1 isa Nothing ? nothing : map(i -> _get_m_expvals(sol[:, i], SaveFuncSSESolve), eachindex(sol))
    m_expvals = _m_expvals isa Nothing ? nothing : stack(_m_expvals, dims = 2)

    expvals =
        _get_expvals(_sol_1, SaveFuncSSESolve) isa Nothing ? nothing :
        dropdims(sum(expvals_all, dims = 2), dims = 2) ./ length(sol)

    return TimeEvolutionStochasticSol(
        ntraj,
        ens_prob.times,
        states,
        expvals,
        expvals, # This is average_expect
        expvals_all,
        m_expvals, # Measurement expectation values
        sol.converged,
        _sol_1.alg,
        _sol_1.prob.kwargs[:abstol],
        _sol_1.prob.kwargs[:reltol],
    )
end
