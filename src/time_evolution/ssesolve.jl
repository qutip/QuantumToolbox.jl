export ssesolveProblem, ssesolveEnsembleProblem, ssesolve

#TODO: Check if works in GPU
function _ssesolve_update_coefficients!(ψ, coefficients, sc_ops)
    _get_en = op -> real(dot(ψ, op, ψ)) #this is en/2: <Sn + Sn'>/2 = Re<Sn>
    @. coefficients[2:end-1] = _get_en(sc_ops) #coefficients of the OperatorSum: Σ Sn * en/2
    coefficients[end] = -sum(x -> x^2, coefficients[2:end-1]) / 2 #this last coefficient is -Σen^2/8
    return nothing
end

function ssesolve_drift!(du, u, p, t)
    _ssesolve_update_coefficients!(u, p.K.coefficients, p.sc_ops)

    mul!(du, p.K, u)

    return nothing
end

function ssesolve_diffusion!(du, u, p, t)
    @inbounds en = @view(p.K.coefficients[2:end-1])

    # du:(H,W). du_reshaped:(H*W,). 
    # H:Hilbert space dimension, W: number of sc_ops
    du_reshaped = reshape(du, :)
    mul!(du_reshaped, p.D, u) #du[:,i] = D[i] * u

    du .-= u .* reshape(en, 1, :) #du[:,i] -= en[i] * u

    return nothing
end

function _ssesolve_prob_func(prob, i, repeat)
    internal_params = prob.p

    global_rng = internal_params.global_rng
    seed = internal_params.seeds[i]
    traj_rng = typeof(global_rng)()
    seed!(traj_rng, seed)

    noise = RealWienerProcess(
        prob.tspan[1],
        zeros(length(internal_params.sc_ops)),
        zeros(length(internal_params.sc_ops)),
        save_everystep = false,
        rng = traj_rng,
    )

    noise_rate_prototype = similar(prob.u0, length(prob.u0), length(internal_params.sc_ops))

    prm = merge(
        internal_params,
        (
            expvals = similar(internal_params.expvals),
            progr = ProgressBar(size(internal_params.expvals, 2), enable = false),
        ),
    )

    return remake(prob, p = prm, noise = noise, noise_rate_prototype = noise_rate_prototype, seed = seed)
end

# Standard output function
_ssesolve_output_func(sol, i) = (sol, false)

# Output function with progress bar update
function _ssesolve_output_func_progress(sol, i)
    next!(sol.prob.p.progr)
    return _ssesolve_output_func(sol, i)
end

# Output function with distributed channel update for progress bar
function _ssesolve_output_func_distributed(sol, i)
    put!(sol.prob.p.progr_channel, true)
    return _ssesolve_output_func(sol, i)
end

_ssesolve_dispatch_output_func() = _ssesolve_output_func
_ssesolve_dispatch_output_func(::ET) where {ET<:Union{EnsembleSerial,EnsembleThreads}} = _ssesolve_output_func_progress
_ssesolve_dispatch_output_func(::EnsembleDistributed) = _ssesolve_output_func_distributed

function _ssesolve_generate_statistics!(sol, i, states, expvals_all)
    sol_i = sol[:, i]
    !isempty(sol_i.prob.kwargs[:saveat]) ?
    states[i] = [QuantumObject(sol_i.u[i], dims = sol_i.prob.p.Hdims) for i in 1:length(sol_i.u)] : nothing

    copyto!(view(expvals_all, i, :, :), sol_i.prob.p.expvals)
    return nothing
end

@doc raw"""
    ssesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector;
        sc_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::StochasticDiffEqAlgorithm=SRA1()
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        rng::AbstractRNG=default_rng(),
        kwargs...)

Generates the SDEProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
    ```math
    d|\psi(t)\rangle = -i K |\psi(t)\rangle dt + \sum_n M_n |\psi(t)\rangle dW_n(t)
    ```

where 
    
```math
    K = \hat{H} + i \sum_n \left(\frac{e_j} C_n - \frac{1}{2} \sum_{j} C_n^\dagger C_n - \frac{e_j^2}{8}\right),
    ```
    ```math
    M_n = C_n - \frac{e_n}{2},
    ```
    and
    ```math
    e_n = \langle C_n + C_n^\dagger \rangle.
    ```

Above, `C_n` is the `n`-th collapse operator and  `dW_j(t)` is the real Wiener increment associated to `C_n`.

# Arguments

- `H::QuantumObject`: The Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: The initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: The time list of the evolution.
- `sc_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: List of stochastic collapse operators ``\{\hat{C}_n\}_n``.
- `alg::StochasticDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::NamedTuple`: The parameters of the system.
- `rng::AbstractRNG`: The random number generator for reproducibility.
- `kwargs...`: The keyword arguments passed to the `SDEProblem` constructor.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (SDE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The `SDEProblem` for the Stochastic Schrödinger time evolution of the system.
"""
function ssesolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    kwargs...,
) where {MT1<:AbstractMatrix,T2}
    H.dims != ψ0.dims && throw(DimensionMismatch("The two quantum objects don't have the same Hilbert dimension."))

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    sc_ops isa Nothing &&
        throw(ArgumentError("The list of collapse operators must be provided. Use sesolveProblem instead."))

    !(H_t isa Nothing) && throw(ArgumentError("Time-dependent Hamiltonians are not currently supported in ssesolve."))

    t_l = convert(Vector{Float64}, tlist) # Convert it into Float64 to avoid type instabilities for StochasticDiffEq.jl

    ϕ0 = get_data(ψ0)

    H_eff = get_data(H - T2(0.5im) * mapreduce(op -> op' * op, +, sc_ops))
    sc_ops2 = get_data.(sc_ops)

    coefficients = [1.0, fill(0.0, length(sc_ops) + 1)...]
    operators = [-1im * H_eff, sc_ops2..., MT1(I(prod(H.dims)))]
    K = OperatorSum(coefficients, operators)
    _ssesolve_update_coefficients!(ϕ0, K.coefficients, sc_ops2)

    D = reduce(vcat, sc_ops2)

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(t_l))
        e_ops2 = MT1[]
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        e_ops2 = get_data.(e_ops)
        is_empty_e_ops = isempty(e_ops)
    end

    p = (
        K = K,
        D = D,
        e_ops = e_ops2,
        sc_ops = sc_ops2,
        expvals = expvals,
        Hdims = H.dims,
        H_t = H_t,
        times = t_l,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = is_empty_e_ops ? t_l : [t_l[end]]
    default_values = (DEFAULT_SDE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, Val(false), t_l, kwargs2)

    tspan = (t_l[1], t_l[end])
    noise = RealWienerProcess(t_l[1], zeros(length(sc_ops)), zeros(length(sc_ops)), save_everystep = false, rng = rng)
    noise_rate_prototype = similar(ϕ0, length(ϕ0), length(sc_ops))
    return SDEProblem{true}(
        ssesolve_drift!,
        ssesolve_diffusion!,
        ϕ0,
        tspan,
        p;
        noise_rate_prototype = noise_rate_prototype,
        noise = noise,
        kwargs3...,
    )
end

@doc raw"""
    ssesolveEnsembleProblem(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector;
        sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::StochasticDiffEqAlgorithm=SRA1()
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        rng::AbstractRNG=default_rng(),
        ntraj::Int=1,
        ensemble_method=EnsembleThreads(),
        prob_func::Function=_mcsolve_prob_func,
        output_func::Function=_ssesolve_dispatch_output_func(ensemble_method),
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Generates the SDE EnsembleProblem for the Stochastic Schrödinger time evolution of a quantum system. This is defined by the following stochastic differential equation:
    
    ```math
    d|\psi(t)\rangle = -i K |\psi(t)\rangle dt + \sum_n M_n |\psi(t)\rangle dW_n(t)
    ```

where 
    
```math
    K = \hat{H} + i \sum_n \left(\frac{e_j} C_n - \frac{1}{2} \sum_{j} C_n^\dagger C_n - \frac{e_j^2}{8}\right),
    ```
    ```math
    M_n = C_n - \frac{e_n}{2},
    ```
    and
    ```math
    e_n = \langle C_n + C_n^\dagger \rangle.
    ```

Above, `C_n` is the `n`-th collapse operator and  `dW_j(t)` is the real Wiener increment associated to `C_n`.

# Arguments

- `H::QuantumObject`: The Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: The initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: The time list of the evolution.
- `sc_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: List of stochastic collapse operators ``\{\hat{C}_n\}_n``.
- `alg::StochasticDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::NamedTuple`: The parameters of the system.
- `rng::AbstractRNG`: The random number generator for reproducibility.
- `ntraj::Int`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use.
- `prob_func::Function`: Function to use for generating the SDEProblem.
- `output_func::Function`: Function to use for generating the output of a single trajectory.
- `progress_bar::Union{Val,Bool}`: Whether to show a progress bar.
- `kwargs...`: The keyword arguments passed to the `SDEProblem` constructor.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (SDE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::EnsembleProblem with SDEProblem`: The Ensemble SDEProblem for the Stochastic Shrödinger time evolution.
"""
function ssesolveEnsembleProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2}
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
            alg = alg,
            e_ops = e_ops,
            H_t = H_t,
            params = merge(params, (global_rng = rng, seeds = seeds)),
            rng = rng,
            kwargs...,
        )

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
    ssesolve(H::QuantumObject,
        ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
        tlist::AbstractVector,
        sc_ops::Union{Nothing, AbstractVector}=nothing;
        alg::StochasticDiffEqAlgorithm=SRA1(),
        e_ops::Union{Nothing,AbstractVector,Tuple}=nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        rng::AbstractRNG=default_rng(),
        ntraj::Int=1,
        ensemble_method=EnsembleThreads(),
        prob_func::Function=_ssesolve_prob_func,
        output_func::Function=_ssesolve_dispatch_output_func(ensemble_method),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...)


Stochastic Schrödinger equation evolution of a quantum system given the system Hamiltonian ``\hat{H}`` and a list of stochadtic collapse (jump) operators ``\{\hat{C}_n\}_n``.
The stochastic evolution of the state ``|\psi(t)\rangle`` is defined by:
    
    ```math
    d|\psi(t)\rangle = -i K |\psi(t)\rangle dt + \sum_n M_n |\psi(t)\rangle dW_n(t)
    ```

where 
    
```math
    K = \hat{H} + i \sum_n \left(\frac{e_j} C_n - \frac{1}{2} \sum_{j} C_n^\dagger C_n - \frac{e_j^2}{8}\right),
    ```
    ```math
    M_n = C_n - \frac{e_n}{2},
    ```
    and
    ```math
    e_n = \langle C_n + C_n^\dagger \rangle.
    ```

Above, `C_n` is the `n`-th collapse operator and  `dW_j(t)` is the real Wiener increment associated to `C_n`.


# Arguments

- `H::QuantumObject`: Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: List of times at which to save the state of the system.
- `sc_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: List of stochastic collapse operators ``\{\hat{C}_n\}_n``.
- `alg::StochasticDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `rng::AbstractRNG`: Random number generator for reproducibility.
- `ntraj::Int`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use.
- `prob_func::Function`: Function to use for generating the SDEProblem.
- `output_func::Function`: Function to use for generating the output of a single trajectory.
- `progress_bar::Union{Val,Bool}`: Whether to show a progress bar.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- `ensemble_method` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`
- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-2` and `abstol=1e-2`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (SDE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionSSESol`: The solution of the time evolution. See also [`TimeEvolutionSSESol`](@ref)
"""
function ssesolve(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    sc_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::StochasticDiffEqAlgorithm = SRA1(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    rng::AbstractRNG = default_rng(),
    ntraj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_dispatch_output_func(ensemble_method),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2}
    progr = ProgressBar(ntraj, enable = getVal(progress_bar))
    progr_channel::RemoteChannel{Channel{Bool}} = RemoteChannel(() -> Channel{Bool}(1))
    @async while take!(progr_channel)
        next!(progr)
    end

    ens_prob = ssesolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        sc_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
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

        expvals_all = Array{ComplexF64}(undef, length(sol), size(_sol_1.prob.p.expvals)...)
        states =
            isempty(_sol_1.prob.kwargs[:saveat]) ? fill(QuantumObject[], length(sol)) :
            Vector{Vector{QuantumObject}}(undef, length(sol))

        foreach(i -> _ssesolve_generate_statistics!(sol, i, states, expvals_all), eachindex(sol))
        expvals = dropdims(sum(expvals_all, dims = 1), dims = 1) ./ length(sol)

        return TimeEvolutionSSESol(
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
