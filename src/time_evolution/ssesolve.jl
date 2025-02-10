export ssesolveProblem, ssesolveEnsembleProblem, ssesolve

# TODO: Merge this with _stochastic_prob_func
function _ssesolve_prob_func(prob, i, repeat)
    internal_params = prob.p

    global_rng = internal_params.global_rng
    seed = internal_params.seeds[i]
    traj_rng = typeof(global_rng)()
    seed!(traj_rng, seed)

    noise = RealWienerProcess!(
        prob.tspan[1],
        zeros(internal_params.n_sc_ops),
        zeros(internal_params.n_sc_ops),
        save_everystep = false,
        rng = traj_rng,
    )

    return remake(prob, noise = noise, seed = seed)
end

#=
    struct DiffusionOperator

A struct to represent the diffusion operator. This is used to perform the diffusion process on N different Wiener processes.
=#
struct DiffusionOperator{T,OpType<:Tuple{Vararg{AbstractSciMLOperator}}} <: AbstractSciMLOperator{T}
    ops::OpType
    function DiffusionOperator(ops::OpType) where {OpType}
        T = mapreduce(eltype, promote_type, ops)
        return new{T,OpType}(ops)
    end
end

@generated function update_coefficients!(L::DiffusionOperator, u, p, t)
    ops_types = L.parameters[2].parameters
    N = length(ops_types)
    quote
        Base.@nexprs $N i -> begin
            update_coefficients!(L.ops[i], u, p, t)
        end

        nothing
    end
end

@generated function LinearAlgebra.mul!(v::AbstractVecOrMat, L::DiffusionOperator, u::AbstractVecOrMat)
    ops_types = L.parameters[2].parameters
    N = length(ops_types)
    quote
        M = length(u)
        S = size(v)
        (S[1] == M && S[2] == $N) || throw(DimensionMismatch("The size of the output vector is incorrect."))
        Base.@nexprs $N i -> begin
            mul!(@view(v[:, i]), L.ops[i], u)
        end
        return v
    end
end

# TODO: Implement the three-argument dot function for SciMLOperators.jl
# Currently, we are assuming a time-independent MatrixOperator
function _ssesolve_update_coeff(u, p, t, op)
    normalize!(u)
    return real(dot(u, op.A, u)) #this is en/2: <Sn + Sn'>/2 = Re<Sn>
end

# Output function with progress bar update
function _ssesolve_output_func_progress(sol, i)
    next!(sol.prob.p.progr)
    return _stochastic_output_func(sol, i)
end

# Output function with distributed channel update for progress bar
function _ssesolve_output_func_distributed(sol, i)
    put!(sol.prob.p.progr_channel, true)
    return _stochastic_output_func(sol, i)
end

_ssesolve_dispatch_output_func(::ET) where {ET<:Union{EnsembleSerial,EnsembleThreads}} = _ssesolve_output_func_progress
_ssesolve_dispatch_output_func(::EnsembleDistributed) = _ssesolve_output_func_distributed

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
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the SDEProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{C}_n - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{C}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{C}_n + \hat{C}_n^\dagger \rangle.
```

Above, ``\hat{C}_n`` is the `n`-th collapse operator and ``dW_n(t)`` is the real Wiener increment associated to ``\hat{C}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
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

- `prob`: The `SDEProblem` for the Stochastic Schrödinger time evolution of the system.
"""
function ssesolveProblem(
    H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
)
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    sc_ops isa Nothing &&
        throw(ArgumentError("The list of collapse operators must be provided. Use sesolveProblem instead."))

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

    p = (progr = progr, times = tlist, Hdims = dims, n_sc_ops = length(sc_ops), params...)

    is_empty_e_ops = (e_ops isa Nothing) ? true : isempty(e_ops)

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (DEFAULT_SDE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_se_me_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2, SaveFuncSSESolve)
    kwargs4 = _ssesolve_add_normalize_cb(kwargs3)

    tspan = (tlist[1], tlist[end])
    noise =
        RealWienerProcess!(tlist[1], zeros(length(sc_ops)), zeros(length(sc_ops)), save_everystep = false, rng = rng)
    noise_rate_prototype = similar(ψ0, length(ψ0), length(sc_ops))
    return SDEProblem{true}(K, D, ψ0, tspan, p; noise_rate_prototype = noise_rate_prototype, noise = noise, kwargs4...)
end

@doc raw"""
    ssesolveEnsembleProblem(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Function = _ssesolve_prob_func,
        output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the SDE EnsembleProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{C}_n - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{C}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{C}_n + \hat{C}_n^\dagger \rangle.
```

Above, ``\hat{C}_n`` is the `n`-th collapse operator and  ``dW_n(t)`` is the real Wiener increment associated to ``\hat{C}_n``. See [Wiseman2009Quantum](@cite) for more details.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `jump_callback`: The Jump Callback type: Discrete or Continuous. The default is `ContinuousLindbladJumpCallback()`, which is more precise.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: Function to use for generating the output of a single trajectory.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
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
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
)
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
        prob_sse = ssesolveProblem(
            H,
            ψ0,
            tlist,
            sc_ops;
            e_ops = e_ops,
            params = merge(params, (global_rng = rng, seeds = seeds)),
            rng = rng,
            progress_bar = Val(false),
            kwargs...,
        )

        # safetycopy is set to true because the K and D functions cannot be currently deepcopied.
        # the memory overhead shouldn't be too large, compared to the safetycopy=false case.
        ensemble_prob = EnsembleProblem(prob_sse, prob_func = prob_func, output_func = output_func, safetycopy = true)

        return ensemble_prob
    catch e
        if ensemble_method isa EnsembleDistributed
            put!(progr_channel, false)
        end
        rethrow(e)
    end
end

@doc raw"""
    ssesolve(
        H::Union{AbstractQuantumObject{OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::StochasticDiffEqAlgorithm = SRA1(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        rng::AbstractRNG = default_rng(),
        ntraj::Int = 1,
        ensemble_method = EnsembleThreads(),
        prob_func::Function = _ssesolve_prob_func,
        output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )


Stochastic Schrödinger equation evolution of a quantum system given the system Hamiltonian ``\hat{H}`` and a list of stochadtic collapse (jump) operators ``\{\hat{C}_n\}_n``.
The stochastic evolution of the state ``|\psi(t)\rangle`` is defined by:
    
```math
d|\psi(t)\rangle = -i \hat{K} |\psi(t)\rangle dt + \sum_n \hat{M}_n |\psi(t)\rangle dW_n(t)
```

where 
    
```math
\hat{K} = \hat{H} + i \sum_n \left(\frac{e_n}{2} \hat{C}_n - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n - \frac{e_n^2}{8}\right),
```
```math
\hat{M}_n = \hat{C}_n - \frac{e_n}{2},
```
and
```math
e_n = \langle \hat{C}_n + \hat{C}_n^\dagger \rangle.
```

Above, ``\hat{C}_n`` is the `n`-th collapse operator and ``dW_n(t)`` is the real Wiener increment associated to ``\hat{C}_n``. See [Wiseman2009Quantum](@cite) for more details.


# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `sc_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm to use for the stochastic differential equation. Default is `SRA1()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` of parameters to pass to the solver.
- `rng`: Random number generator for reproducibility.
- `ntraj`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use. Default to `EnsembleThreads()`.
- `prob_func`: Function to use for generating the ODEProblem.
- `output_func`: Function to use for generating the output of a single trajectory.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
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
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
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
        kwargs...,
    )

    return ssesolve(ens_prob; alg = alg, ntraj = ntraj, ensemble_method = ensemble_method)
end

function ssesolve(
    ens_prob::EnsembleProblem;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
)
    # Stop the async task if an error occurs
    try
        sol = solve(ens_prob, alg, ensemble_method, trajectories = ntraj)

        if ensemble_method isa EnsembleDistributed
            put!(sol[:, 1].prob.p.progr_channel, false)
        end

        _sol_1 = sol[:, 1]
        _expvals_sol_1 = _se_me_sse_get_expvals(_sol_1)

        normalize_states = Val(false)
        dims = _sol_1.prob.p.Hdims
        _expvals_all =
            _expvals_sol_1 isa Nothing ? nothing : map(i -> _se_me_sse_get_expvals(sol[:, i]), eachindex(sol))
        expvals_all = _expvals_all isa Nothing ? nothing : stack(_expvals_all)
        states = map(i -> _normalize_state!.(sol[:, i].u, Ref(dims), normalize_states), eachindex(sol))

        expvals =
            _se_me_sse_get_expvals(_sol_1) isa Nothing ? nothing :
            dropdims(sum(expvals_all, dims = 3), dims = 3) ./ length(sol)

        return TimeEvolutionStochasticSol(
            ntraj,
            _sol_1.prob.p.times,
            states,
            expvals,
            expvals_all,
            sol.converged,
            _sol_1.alg,
            _sol_1.prob.kwargs[:abstol],
            _sol_1.prob.kwargs[:reltol],
        )
    catch e
        if ensemble_method isa EnsembleDistributed
            put!(ens_prob.prob.p.progr_channel, false)
        end
        rethrow(e)
    end
end
