export ssesolveProblem, ssesolveEnsembleProblem, ssesolve

#TODO: Check if works in GPU
function _ssesolve_update_coefficients!(ψ, coefficients, c_ops)
    _get_en = op -> real(dot(ψ, op, ψ)) #this is en/2: <Sn + Sn'>/2 = Re<Sn>
    @. coefficients[2:end-1] = _get_en(c_ops) #coefficients of the OperatorSum: Σ Sn * en/2
    coefficients[end] = -sum(x -> x^2, coefficients[2:end-1]) / 2 #this last coefficient is -Σen^2/8
    return nothing
end

function ssesolve_drift!(du, u, p, t)
    _ssesolve_update_coefficients!(u, p.K.coefficients, p.c_ops)

    mul!(du, p.K, u)

    return nothing
end

function ssesolve_diffusion!(du, u, p, t)
    @inbounds en = @view(p.K.coefficients[2:end-1])

    # du:(H,W). du_reshaped:(H*W,). 
    # H:Hilbert space dimension, W: number of c_ops
    du_reshaped = reshape(du, :)
    mul!(du_reshaped, p.D, u) #du[:,i] = D[i] * u

    du .-= u .* reshape(en, 1, :) #du[:,i] -= en[i] * u

    return nothing
end

function _ssesolve_prob_func(prob, i, repeat)
    internal_params = prob.p

    noise_rate_prototype = similar(prob.u0, length(prob.u0), length(internal_params.c_ops))
    noise = RealWienerProcess!(prob.tspan[1], zeros(length(internal_params.c_ops)))

    prm = merge(
        internal_params,
        (
            expvals = similar(internal_params.expvals),
            progr = ProgressBar(size(internal_params.expvals, 2), enable = false),
        ),
    )

    return remake(prob, p = prm, noise_rate_prototype = noise_rate_prototype, noise = noise)
end

_ssesolve_output_func(sol, i) = (sol, false)

function _ssesolve_generate_statistics!(sol, i, states, expvals_all)
    sol_i = sol[:, i]
    !isempty(sol_i.prob.kwargs[:saveat]) ?
    states[i] = [QuantumObject(sol_i.u[i], dims = sol_i.prob.p.Hdims) for i in 1:length(sol_i.u)] : nothing

    copyto!(view(expvals_all, i, :, :), sol_i.prob.p.expvals)
    return nothing
end

function ssesolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::StochasticDiffEq.StochasticDiffEqAlgorithm = EM(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix}
    H.dims != ψ0.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    !(H_t isa Nothing) && throw(ArgumentError("Time-dependent Hamiltonians are not currently supported in ssesolve."))

    progress_bar_val = makeVal(progress_bar)

    t_l = convert(Vector{Float64}, tlist) # Convert it into Float64 to avoid type instabilities for OrdinaryDiffEq.jl

    ϕ0 = get_data(ψ0)

    H_eff = get_data(H - T2(0.5im) * mapreduce(op -> op' * op, +, c_ops))
    c_ops2 = get_data.(c_ops)

    coefficients = [1.0, fill(0.0, length(c_ops) + 1)...]
    operators = [-1im * H_eff, c_ops2..., I(prod(H.dims))]
    K = OperatorSum(coefficients, operators)
    _ssesolve_update_coefficients!(ϕ0, K.coefficients, c_ops2)

    D = vcat(c_ops2...)

    progr = ProgressBar(length(t_l), enable = getVal(progress_bar_val))

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
        c_ops = c_ops2,
        expvals = expvals,
        progr = progr,
        Hdims = H.dims,
        H_t = H_t,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = e_ops isa Nothing ? t_l : [t_l[end]]
    default_values = (DEFAULT_SDE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, progress_bar_val, t_l, kwargs2)

    tspan = (t_l[1], t_l[end])
    noise = RealWienerProcess!(t_l[1], zeros(length(c_ops)))
    noise_rate_prototype = similar(ϕ0, length(ϕ0), length(c_ops))
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

function ssesolveEnsembleProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::StochasticDiffEq.StochasticDiffEqAlgorithm = EM(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_output_func,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix}
    prob_sse = ssesolveProblem(H, ψ0, tlist, c_ops; alg = alg, e_ops = e_ops, H_t = H_t, params = params, kwargs...)

    ensemble_prob = EnsembleProblem(prob_sse, prob_func = prob_func, output_func = output_func, safetycopy = false)

    return ensemble_prob
end

function ssesolve(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    tlist::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[];
    alg::StochasticDiffEq.StochasticDiffEqAlgorithm = EM(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    n_traj::Int = 1,
    ensemble_method = EnsembleThreads(),
    prob_func::Function = _ssesolve_prob_func,
    output_func::Function = _ssesolve_output_func,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,Tc<:AbstractMatrix}
    ens_prob = ssesolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        prob_func = prob_func,
        output_func = output_func,
        kwargs...,
    )

    return ssesolve(ens_prob; alg = alg, n_traj = n_traj, ensemble_method = ensemble_method)
end

function ssesolve(
    ens_prob::EnsembleProblem;
    alg::StochasticDiffEq.StochasticDiffEqAlgorithm = EM(),
    n_traj::Int = 1,
    ensemble_method = EnsembleThreads(),
)
    sol = solve(ens_prob, alg, ensemble_method, trajectories = n_traj)
    _sol_1 = sol[:, 1]

    expvals_all = Array{ComplexF64}(undef, length(sol), size(_sol_1.prob.p.expvals)...)
    states =
        isempty(_sol_1.prob.kwargs[:saveat]) ? fill(QuantumObject[], length(sol)) :
        Vector{Vector{QuantumObject}}(undef, length(sol))

    foreach(i -> _ssesolve_generate_statistics!(sol, i, states, expvals_all), eachindex(sol))
    expvals = dropdims(sum(expvals_all, dims = 1), dims = 1) ./ length(sol)

    return TimeEvolutionSSESol(
        n_traj,
        _sol_1.t,
        states,
        expvals,
        expvals_all,
        sol.converged,
        _sol_1.alg,
        _sol_1.prob.kwargs[:abstol],
        _sol_1.prob.kwargs[:reltol],
    )
end
