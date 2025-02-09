export smesolveProblem, smesolveEnsembleProblem, smesolve

_smesolve_generate_state(u, dims) = QuantumObject(vec2mat(u), type = Operator, dims = dims)

function _smesolve_update_coeff(u, p, t, op_vec)
    return real(dot(u, op_vec)) / 2 #this is Tr[Sn * ρ + ρ * Sn']
end

_smesolve_ScalarOperator(op_vec) =
    ScalarOperator(one(eltype(op_vec)), (a, u, p, t) -> -_smesolve_update_coeff(u, p, t, op_vec))

@doc raw"""
    smesolveProblem(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d| \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_n \mathcal{D}[\hat{C}_n] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,

Above, ``\hat{C}_n`` represent the operators related to pure dissipation, while ``\hat{S}_n`` are the measurement operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref) or a [`Operator`](@ref).
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of measurement collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the `SDEProblem` for the Stochastic Master Equation time evolution.
"""
function smesolveProblem(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
)
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    isnothing(sc_ops) &&
        throw(ArgumentError("The list of measurement collapse operators must be provided. Use mesolveProblem instead."))

    tlist = _check_tlist(tlist, _FType(ψ0))

    L_evo = _mesolve_make_L_QobjEvo(H, c_ops) + _mesolve_make_L_QobjEvo(nothing, sc_ops)
    check_dimensions(L_evo, ψ0)
    dims = L_evo.dimensions

    T = Base.promote_eltype(L_evo, ψ0)
    ρ0 = sparse_to_dense(_CType(T), mat2vec(ket2dm(ψ0).data)) # Convert it to dense vector with complex element type

    progr = ProgressBar(length(tlist), enable = getVal(progress_bar))

    sc_ops_evo_data = Tuple(map(get_data ∘ QobjEvo, sc_ops))

    K = get_data(L_evo)

    Id = I(prod(dims))
    D_l = map(sc_ops_evo_data) do op
        # TODO: Implement the three-argument dot function for SciMLOperators.jl
        # Currently, we are assuming a time-independent MatrixOperator
        op_vec = mat2vec(adjoint(op.A))
        return _spre(op, Id) + _spost(op', Id) + _smesolve_ScalarOperator(op_vec) * IdentityOperator(prod(dims)^2)
    end
    D = DiffusionOperator(D_l)

    p = (progr = progr, times = tlist, Hdims = dims, n_sc_ops = length(sc_ops), params...)

    is_empty_e_ops = (e_ops isa Nothing) ? true : isempty(e_ops)

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (DEFAULT_SDE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_se_me_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2, SaveFuncMESolve)

    tspan = (tlist[1], tlist[end])
    noise =
        RealWienerProcess!(tlist[1], zeros(length(sc_ops)), zeros(length(sc_ops)), save_everystep = false, rng = rng)
    noise_rate_prototype = similar(ρ0, length(ρ0), length(sc_ops))
    prob = SDEProblem{true}(K, D, ρ0, tspan, p; noise_rate_prototype = noise_rate_prototype, noise = noise, kwargs3...)

    return TimeEvolutionProblem(prob, tlist, dims)
end

@doc raw"""
    smesolveEnsembleProblem(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d| \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_n \mathcal{D}[\hat{C}_n] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,

Above, ``\hat{C}_n`` represent the operators related to pure dissipation, while ``\hat{S}_n`` are the measurement operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref) or a [`Operator`](@ref).
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of measurement collapse operators ``\{\hat{S}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the Ensemble `SDEProblem` for the Stochastic Master Equation time evolution.
"""
function smesolveEnsembleProblem(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
)
    _prob_func =
        isnothing(prob_func) ? _ensemble_dispatch_prob_func(rng, ntraj, tlist, _stochastic_prob_func) : prob_func
    _output_func =
        output_func isa Nothing ?
        _ensemble_dispatch_output_func(ensemble_method, progress_bar, ntraj, _stochastic_output_func) : output_func

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
    smesolve(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::StochasticDiffEqAlgorithm = SRA1(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Union{Function, Nothing} = nothing,
        output_func::Union{Tuple,Nothing} = nothing,
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Stochastic Master Equation time evolution of an open quantum system. This is defined by the following stochastic differential equation:
    
```math
d| \rho (t) = -i [\hat{H}, \rho(t)] dt + \sum_n \mathcal{D}[\hat{C}_n] \rho(t) dt + \sum_n \mathcal{D}[\hat{S}_n] \rho(t) dt + \sum_n \mathcal{H}[\hat{S}_n] \rho(t) dW_n(t),
```

where

```math
\mathcal{D}[\hat{O}] \rho = \hat{O} \rho \hat{O}^\dagger - \frac{1}{2} \{\hat{O}^\dagger \hat{O}, \rho\},
```

is the Lindblad superoperator, and

```math
\mathcal{H}[\hat{O}] \rho = \hat{O} \rho + \rho \hat{O}^\dagger - \mathrm{Tr}[\hat{O} \rho + \rho \hat{O}^\dagger] \rho,

Above, ``\hat{C}_n`` represent the operators related to pure dissipation, while ``\hat{S}_n`` are the measurement operators. The ``dW_n(t)`` term is the real Wiener increment associated to ``\hat{S}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref) or a [`Operator`](@ref).
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `sc_ops`: List of measurement collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm to use for the stochastic differential equation. Default is `SRA1()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `ProgressBar` object, and the (optional) `RemoteChannel` object.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionSMESol`: The solution of the time evolution. See also [`TimeEvolutionSMESol`](@ref).
"""
function smesolve(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
)
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
        ensemble_method = ensemble_method,
        prob_func = prob_func,
        output_func = output_func,
        progress_bar = progress_bar,
        kwargs...,
    )

    return smesolve(ensemble_prob, alg, ntraj, ensemble_method)
end

function smesolve(
    ens_prob::TimeEvolutionProblem,
    alg::StochasticDiffEqAlgorithm = SRA1(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
)
    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemble_method, ntraj)

    _sol_1 = sol[:, 1]
    _expvals_sol_1 = _se_me_sse_get_expvals(_sol_1)

    normalize_states = Val(false)
    dims = ens_prob.dimensions
    _expvals_all = _expvals_sol_1 isa Nothing ? nothing : map(i -> _se_me_sse_get_expvals(sol[:, i]), eachindex(sol))
    expvals_all = _expvals_all isa Nothing ? nothing : stack(_expvals_all)
    states = map(i -> _smesolve_generate_state.(sol[:, i].u, Ref(dims)), eachindex(sol))

    expvals =
        _se_me_sse_get_expvals(_sol_1) isa Nothing ? nothing :
        dropdims(sum(expvals_all, dims = 3), dims = 3) ./ length(sol)

    return TimeEvolutionSMESol(
        ntraj,
        ens_prob.times,
        states,
        expvals,
        expvals_all,
        sol.converged,
        _sol_1.alg,
        _sol_1.prob.kwargs[:abstol],
        _sol_1.prob.kwargs[:reltol],
    )
end
