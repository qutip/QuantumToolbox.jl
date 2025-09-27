export smesolveProblem, smesolveEnsembleProblem, smesolve

_smesolve_generate_state(u, dims, isoperket::Val{false}) = QuantumObject(vec2mat(u), type = Operator(), dims = dims)
_smesolve_generate_state(u, dims, isoperket::Val{true}) = QuantumObject(u, type = OperatorKet(), dims = dims)

function _smesolve_update_coeff(u, p, t, op_vec)
    return 2 * real(dot(op_vec, u)) #this is Tr[Sn * ρ + ρ * Sn']
end

_smesolve_ScalarOperator(op_vec) =
    ScalarOperator(one(eltype(op_vec)), (a, u, p, t) -> -_smesolve_update_coeff(u, p, t, op_vec))

@doc raw"""
    smesolveProblem(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject,
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        progress_bar::Union{Val,Bool} = Val(true),
        store_measurement::Union{Val, Bool} = Val(false),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_i \mathcal{D}[\hat{C}_i] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,
```

Above, ``\hat{C}_i`` represent the collapse operators related to pure dissipation, while ``\hat{S}_n`` are the stochastic collapse operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_i\}_i``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector`, a `Tuple` or a [`AbstractQuantumObject`](@ref). It is recommended to use the last case when only one operator is provided.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `store_measurement`: Whether to store the measurement expectation values. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

!!! tip "Performance Tip"
    When `sc_ops` contains only a single operator, it is recommended to pass only that operator as the argument. This ensures that the stochastic noise is diagonal, making the simulation faster.

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the `SDEProblem` for the Stochastic Master Equation time evolution.
"""
function smesolveProblem(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    progress_bar::Union{Val,Bool} = Val(true),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
) where {StateOpType<:Union{Ket,Operator,OperatorKet}}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    isnothing(sc_ops) &&
        throw(ArgumentError("The list of stochastic collapse operators must be provided. Use mesolveProblem instead."))
    sc_ops_list = _make_c_ops_list(sc_ops) # If it is an AbstractQuantumObject but we need to iterate
    sc_ops_isa_Qobj = sc_ops isa AbstractQuantumObject # We can avoid using non-diagonal noise if sc_ops is just an AbstractQuantumObject

    tlist = _check_tlist(tlist, _float_type(ψ0))

    L_evo = _mesolve_make_L_QobjEvo(H, c_ops) + _mesolve_make_L_QobjEvo(nothing, sc_ops_list)
    check_dimensions(L_evo, ψ0)
    dims = L_evo.dimensions

    T = Base.promote_eltype(L_evo, ψ0)
    ρ0 = if isoperket(ψ0) # Convert it to dense vector with complex element type
        to_dense(_complex_float_type(T), copy(ψ0.data))
    else
        to_dense(_complex_float_type(T), mat2vec(ket2dm(ψ0).data))
    end

    progr = ProgressBar(length(tlist), enable = getVal(progress_bar))

    sc_ops_evo_data = Tuple(map(get_data ∘ QobjEvo, sc_ops_list))

    K = get_data(L_evo)

    Id = I(prod(dims))
    Id_op = IdentityOperator(prod(dims)^2)
    D_l = map(sc_ops_evo_data) do op
        # TODO: # Currently, we are assuming a time-independent MatrixOperator
        # Also, the u state may become non-hermitian, so Tr[Sn * ρ + ρ * Sn'] != real(Tr[Sn * ρ]) / 2
        op_vec = mat2vec(adjoint(op.A))
        return _spre(op, Id) + _spost(op', Id) + _smesolve_ScalarOperator(op_vec) * Id_op
    end
    D = DiffusionOperator(D_l)

    kwargs2 = _merge_saveat(tlist, e_ops, DEFAULT_SDE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_stochastic_kwargs(
        e_ops,
        sc_ops_list,
        makeVal(progress_bar),
        tlist,
        makeVal(store_measurement),
        kwargs2,
        SaveFuncSMESolve,
    )

    tspan = (tlist[1], tlist[end])
    noise = _make_noise(tspan[1], sc_ops, makeVal(store_measurement), rng)
    noise_rate_prototype = sc_ops_isa_Qobj ? nothing : similar(ρ0, length(ρ0), length(sc_ops_list))
    prob = SDEProblem{true}(
        K,
        D,
        ρ0,
        tspan,
        params;
        noise_rate_prototype = noise_rate_prototype,
        noise = noise,
        kwargs3...,
    )

    return TimeEvolutionProblem(prob, tlist, dims, (isoperket = Val(isoperket(ψ0)),))
end

@doc raw"""
    smesolveEnsembleProblem(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject,
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        store_measurement::Union{Val,Bool} = Val(false),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_i \mathcal{D}[\hat{C}_i] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,
```

Above, ``\hat{C}_i`` represent the collapse operators related to pure dissipation, while ``\hat{S}_n`` are the stochastic collapse operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_i\}_i``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector`, a `Tuple` or a [`AbstractQuantumObject`](@ref). It is recommended to use the last case when only one operator is provided.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use. Default is `500`.
- `ensemblealg`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the SDEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `store_measurement`: Whether to store the measurement expectation values. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

!!! tip "Performance Tip"
    When `sc_ops` contains only a single operator, it is recommended to pass only that operator as the argument. This ensures that the stochastic noise is diagonal, making the simulation faster.

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the Ensemble `SDEProblem` for the Stochastic Master Equation time evolution.
"""
function smesolveEnsembleProblem(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 500,
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
) where {StateOpType<:Union{Ket,Operator,OperatorKet}}
    _prob_func =
        isnothing(prob_func) ?
        _ensemble_dispatch_prob_func(
            rng,
            ntraj,
            tlist,
            _stochastic_prob_func;
            sc_ops = sc_ops,
            store_measurement = makeVal(store_measurement),
        ) : prob_func
    _output_func =
        output_func isa Nothing ?
        _ensemble_dispatch_output_func(ensemblealg, progress_bar, ntraj, _stochastic_output_func) : output_func

    prob_sme = smesolveProblem(
        H,
        ψ0,
        tlist,
        c_ops,
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
        merge(prob_sme.kwargs, (progr = _output_func[2], channel = _output_func[3])),
    )

    return ensemble_prob
end

@doc raw"""
    smesolve(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject,
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
        alg::Union{Nothing,StochasticDiffEqAlgorithm} = nothing,
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        keep_runs_results::Union{Val,Bool} = Val(false),
        store_measurement::Union{Val,Bool} = Val(false),
        kwargs...,
    )

Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_i \mathcal{D}[\hat{C}_i] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,
```

Above, ``\hat{C}_i`` represent the collapse operators related to pure dissipation, while ``\hat{S}_n`` are the stochastic co operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_i\}_i``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of stochastic collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector`, a `Tuple` or a [`AbstractQuantumObject`](@ref). It is recommended to use the last case when only one operator is provided.
- `alg`: The algorithm to use for the stochastic differential equation. Default is `SRIW1()` if `sc_ops` is an [`AbstractQuantumObject`](@ref) (diagonal noise), and `SRA2()` otherwise (non-diagonal noise).
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NullParameters` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use. Default is `500`.
- `ensemblealg`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the SDEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `keep_runs_results`: Whether to save the results of each trajectory. Default to `Val(false)`.
- `store_measurement`: Whether to store the measurement expectation values. Default is `Val(false)`.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

!!! tip "Performance Tip"
    When `sc_ops` contains only a single operator, it is recommended to pass only that operator as the argument. This ensures that the stochastic noise is diagonal, making the simulation faster.

# Returns

- `sol::TimeEvolutionStochasticSol`: The solution of the time evolution. See [`TimeEvolutionStochasticSol`](@ref).
"""
function smesolve(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple,AbstractQuantumObject} = nothing;
    alg::Union{Nothing,StochasticDiffEqAlgorithm} = nothing,
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 500,
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    keep_runs_results::Union{Val,Bool} = Val(false),
    store_measurement::Union{Val,Bool} = Val(false),
    kwargs...,
) where {StateOpType<:Union{Ket,Operator,OperatorKet}}
    ensemble_prob = smesolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        c_ops,
        sc_ops;
        e_ops = e_ops,
        params = params,
        rng = rng,
        ntraj = ntraj,
        ensemblealg = ensemblealg,
        prob_func = prob_func,
        output_func = output_func,
        progress_bar = progress_bar,
        store_measurement = makeVal(store_measurement),
        kwargs...,
    )

    sc_ops_isa_Qobj = sc_ops isa AbstractQuantumObject # We can avoid using non-diagonal noise if sc_ops is just an AbstractQuantumObject

    if isnothing(alg)
        alg = sc_ops_isa_Qobj ? SRIW1() : SRA2()
    end

    return smesolve(ensemble_prob, alg, ntraj, ensemblealg, makeVal(keep_runs_results))
end

function smesolve(
    ens_prob::TimeEvolutionProblem,
    alg::StochasticDiffEqAlgorithm = SRA2(),
    ntraj::Int = 500,
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    keep_runs_results = Val(false),
)
    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemblealg, ntraj)

    _sol_1 = sol[:, 1]
    _expvals_sol_1 = _get_expvals(_sol_1, SaveFuncMESolve)
    _m_expvals_sol_1 = _get_m_expvals(_sol_1, SaveFuncSMESolve)

    dims = ens_prob.dimensions
    _expvals_all =
        _expvals_sol_1 isa Nothing ? nothing : map(i -> _get_expvals(sol[:, i], SaveFuncMESolve), eachindex(sol))
    expvals_all = _expvals_all isa Nothing ? nothing : stack(_expvals_all, dims = 2) # Stack on dimension 2 to align with QuTiP

    # stack to transform Vector{Vector{QuantumObject}} -> Matrix{QuantumObject}
    states_all = stack(
        map(i -> _smesolve_generate_state.(sol[:, i].u, Ref(dims), ens_prob.kwargs.isoperket), eachindex(sol)),
        dims = 1,
    )

    _m_expvals =
        _m_expvals_sol_1 isa Nothing ? nothing : map(i -> _get_m_expvals(sol[:, i], SaveFuncSMESolve), eachindex(sol))
    m_expvals = _m_expvals isa Nothing ? nothing : stack(_m_expvals, dims = 2) # Stack on dimension 2 to align with QuTiP

    kwargs = NamedTuple(_sol_1.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionStochasticSol(
        ntraj,
        ens_prob.times,
        _sol_1.t,
        _store_multitraj_states(states_all, keep_runs_results),
        _store_multitraj_expect(expvals_all, keep_runs_results),
        m_expvals, # Measurement expectation values
        sol.converged,
        _sol_1.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end
