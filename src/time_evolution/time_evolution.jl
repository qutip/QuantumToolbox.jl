abstract type LiouvillianSolver end
struct LiouvillianDirectSolver <: LiouvillianSolver 
    tol::Real
end

abstract type SteadyStateSolver end
abstract type SteadyStateDirectSolver <: SteadyStateSolver end

struct TimeEvolutionSol{TT<:Vector{<:Real}, TS<:AbstractVector, TE<:Matrix{ComplexF64}}
    times::TT
    states::TS
    expect::TE
end

struct TimeEvolutionMCSol{TT<:Vector{<:Vector{<:Real}}, TS<:AbstractVector, TE<:Matrix{ComplexF64}, 
                TEA<:Array{ComplexF64, 3}, TJT<:Vector{<:Vector{<:Real}}, TJW<:Vector{<:Vector{<:Integer}}}
    times::TT
    states::TS
    expect::TE
    expect_all::TEA
    jump_times::TJT
    jump_which::TJW
end

LiouvillianDirectSolver(;tol=1e-8) = LiouvillianDirectSolver(tol)

function _save_func_sesolve(integrator)
    internal_params = integrator.p
    progr = internal_params.progr
    e_ops = internal_params.e_ops
    expvals = internal_params.expvals

    if !isempty(e_ops)
        ψ = integrator.u
        _expect = op -> dot(ψ, op, ψ)
        @. expvals[:, progr.counter+1] = _expect.(e_ops)
    end
    next!(progr)
    u_modified!(integrator, false)
end

function _save_func_mesolve(integrator)
    internal_params = integrator.p
    progr = internal_params.progr
    e_ops = internal_params.e_ops
    expvals = internal_params.expvals

    if !isempty(e_ops)
        # This is equivalent to tr(op * ρ), when both are matrices.
        # The advantage of using this convention is that I don't need
        # to reshape u to make it a matrix, but I reshape the e_ops once.
        
        ρ = integrator.u
        _expect = op -> dot(op, ρ)
        @. expvals[:, progr.counter+1] = _expect.(e_ops)
    end
    next!(progr)
    u_modified!(integrator, false)
end

function _save_func_mcsolve(integrator)
    internal_params = integrator.p
    save_it = internal_params.save_it
    e_ops = internal_params.e_ops
    expvals = internal_params.expvals
    cache_mc = internal_params.cache_mc

    if !isempty(e_ops)
        cache_mc .= integrator.u
        normalize!(cache_mc)
        ψ = cache_mc
        _expect = op -> dot(ψ, op, ψ)
        @. expvals[:, save_it[]+1] = _expect.(e_ops)
    end
    save_it[] += 1
    u_modified!(integrator, false)
end

function LindbladJumpAffect!(integrator)
    internal_params = integrator.p
    c_ops = internal_params.c_ops
    cache_mc = internal_params.cache_mc
    weights_mc = internal_params.weights_mc
    cumsum_weights_mc = internal_params.cumsum_weights_mc
    random_n = internal_params.random_n
    jump_times = internal_params.jump_times
    jump_which = internal_params.jump_which
    ψ = integrator.u

    @inbounds for i in eachindex(weights_mc)
        mul!(cache_mc, c_ops[i], ψ)
        weights_mc[i] = real(dot(cache_mc, cache_mc))
    end
    cumsum!(cumsum_weights_mc, weights_mc)
    collaps_idx = getindex(1:length(weights_mc), findfirst(>(rand()*sum(weights_mc)), cumsum_weights_mc))
    mul!(cache_mc, c_ops[collaps_idx], ψ)
    normalize!(cache_mc)
    integrator.u .= cache_mc

    push!(jump_times, integrator.t)
    push!(jump_which, collaps_idx)
    random_n[] = rand()
end

function ContinuousLindbladJumpCallback(interp_points::Int=0; affect_func::Function=LindbladJumpAffect!)
    LindbladJumpCondition(u, t, integrator) = integrator.p.random_n[] - real(dot(u, u))

    ContinuousCallback(LindbladJumpCondition, affect_func, nothing, interp_points=interp_points, save_positions=(false, false))
end

function DiscreteLindbladJumpCallback(;affect_func::Function=LindbladJumpAffect!)
    LindbladJumpCondition(u, t, integrator) = real(dot(u, u)) < integrator.p.random_n[]

    DiscreteCallback(LindbladJumpCondition, affect_func, save_positions=(false, false))
end

function _mcsolve_prob_func(prob, i, repeat)
    internal_params = prob.p

    prm = merge(internal_params, (U = deepcopy(internal_params.U), e_ops = deepcopy(internal_params.e_ops), 
                c_ops = deepcopy(internal_params.c_ops), expvals = similar(internal_params.expvals), 
                cache_mc = similar(internal_params.cache_mc), weights_mc = similar(internal_params.weights_mc), 
                cumsum_weights_mc = similar(internal_params.weights_mc), random_n = Ref(rand()), save_it = Ref{Int32}(0),
                jump_times = similar(internal_params.jump_times), jump_which = similar(internal_params.jump_which)))

    remake(prob, p=prm)
end

function _mcsolve_output_func(sol, i)
    internal_params = sol.prob.p
    put!(internal_params.progr_channel, true)
    (sol, false)
end

function _mcsolve_generate_statistics(sol, i, times, states, expvals_all, jump_times, jump_which)
    sol_i = sol[i]
    sol_u = haskey(sol_i.prob.kwargs, :save_idxs) ? sol_i.u : QuantumObject.(sol_i.u, dims=sol_i.prob.p.Hdims)

    expvals_all[i, :, :] .= sol_i.prob.p.expvals
    push!(times, sol_i.t)
    push!(states, sol_u)
    push!(jump_times, sol_i.prob.p.jump_times)
    push!(jump_which, sol_i.prob.p.jump_which)
end





#######################################



"""
    sesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5()
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        kwargs...)

Generates the ODEProblem for the Schrödinger time evolution of a quantum system.

# Arguments
- `H::QuantumObject`: The Hamiltonian of the system.
- `ψ0::QuantumObject`: The initial state of the system.
- `t_l::AbstractVector`: The time list of the evolution.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::AbstractVector`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::Dict{Symbol, Any}`: The parameters of the system.
- `progress::Bool`: Whether to show the progress bar.
- `kwargs...`: The keyword arguments passed to the `ODEProblem` constructor.

# Returns
- `prob`: The `ODEProblem` for the Schrödinger time evolution of the system.
"""
function sesolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    kwargs...) where {T1,T2}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ϕ0 = get_data(ψ0)
    U = -1im * get_data(H)

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = length(e_ops) == 0 ? Vector{Matrix{T1}}([]) : get_data.(e_ops)
    p = merge((U = U, e_ops = e_ops2, expvals = expvals, progr = progr, Hdims = H.dims), params)

    kwargs2 = kwargs
    if !isempty(e_ops) || progress
        cb1 = PresetTimeCallback(t_l, _save_func_sesolve, save_positions=(false, false))
        kwargs2 = merge(kwargs2, haskey(kwargs2, :callback) ? 
                    Dict(:callback => CallbackSet(cb1, kwargs2[:callback])) : Dict(:callback => cb1))
    end
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])

    if isa(alg, OrdinaryDiffEq.OrdinaryDiffEqExponentialAlgorithm)
        is_time_dependent && error("The Liouvillian must be time independent when using LinearExponential algorithm.")
        kwargs2 = merge(kwargs2, Dict(:dt => t_l[2] - t_l[1]))
        A = DiffEqArrayOperator(U)

        return ODEProblem(A, ϕ0, tspan, p; kwargs2...)
    else
        if !is_time_dependent
            dudt! = (du, u, p, t) -> mul!(du, p.U, u)
        else
            dudt! = (du, u, p, t) -> mul!(du, p.U - 1im * H_t(t).data, u)
        end

        return ODEProblem(dudt!, ϕ0, tspan, p; kwargs2...)
    end
end


"""
    mesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector=[];
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        kwargs...)

Generates the ODEProblem for the master equation time evolution of an open quantum system.

# Arguments
- `H::QuantumObject`: The Hamiltonian or the Liouvillian of the system.
- `ψ0::QuantumObject`: The initial state of the system.
- `t_l::AbstractVector`: The time list of the evolution.
- `c_ops::AbstractVector=[]`: The list of the collapse operators.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5()`: The algorithm used for the time evolution.
- `e_ops::AbstractVector=[]`: The list of the operators for which the expectation values are calculated.
- `H_t::Union{Nothing,Function}=nothing`: The time-dependent Hamiltonian or Liouvillian.
- `params::Dict{Symbol, Any}=Dict{Symbol, Any}()`: The parameters of the time evolution.
- `progress::Bool=true`: Whether to show the progress bar.
- `kwargs...`: The keyword arguments for the ODEProblem.

# Returns
- `prob::ODEProblem`: The ODEProblem for the master equation time evolution.
"""
function mesolveProblem(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector=[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    kwargs...) where {T1,T2,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ρ0 = isket(ψ0) ? mat2vec(ket2dm(ψ0).data) : mat2vec(ψ0.data)
    L = liouvillian(H, c_ops).data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = length(e_ops) == 0 ? Vector{Vector{T1}}([]) : map(op -> mat2vec(get_data(op)'), e_ops)
    p = merge((L = L, e_ops = e_ops2, expvals = expvals, progr = progr, Hdims = H.dims), params)

    kwargs2 = kwargs
    if !isempty(e_ops) || progress
        cb1 = PresetTimeCallback(t_l, _save_func_mesolve, save_positions=(false, false))
        kwargs2 = merge(kwargs2, haskey(kwargs2, :callback) ? 
                    Dict(:callback => CallbackSet(cb1, kwargs2[:callback])) : Dict(:callback => cb1))
    end
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])
    if isa(alg, OrdinaryDiffEq.OrdinaryDiffEqExponentialAlgorithm)
        is_time_dependent && error("The Liouvillian must be time independent when using LinearExponential algorithm.")
        kwargs2 = merge(kwargs2, Dict(:dt => t_l[2] - t_l[1]))
        A = DiffEqArrayOperator(L)

        return ODEProblem(A, ρ0, tspan, p; kwargs2...)
    else  
        if !is_time_dependent
            dudt! = (du, u, p, t) -> mul!(du, p.L, u)
        else
            if isoper(H_t(0.0))
                @warn string("To speed up the calculation, it is always better to define ",
                    "the time-dependent part as a SuperOperator, and not as an Operator.") maxlog = 1
                dudt! = (du, u, p, t) -> mul!(du, p.L + liouvillian(H_t(t)).data, u)
            else
                dudt! = (du, u, p, t) -> mul!(du, p.L + H_t(t).data, u)
            end
        end

        return ODEProblem(dudt!, ρ0, tspan, p; kwargs2...)
    end
end


"""
    mcsolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        jump_interp_pts::Integer=-1,
        kwargs...)

Generates the ODEProblem for a single trajectory of the Monte Carlo wave function
time evolution of an open quantum system.

# Arguments
- `H::QuantumObject`: Hamiltonian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::AbstractVector`: List of collapse operators.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function}`: Time-dependent part of the Hamiltonian.
- `params::Dict{Symbol, Any}`: Dictionary of parameters to pass to the solver.
- `jump_interp_pts::Integer`: Number of points to use for interpolation of the jump times.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns
- `prob::ODEProblem`: The ODEProblem for the Monte Carlo wave function time evolution.

# Notes
When `jump_interp_pts` is set to -1, a `DiscreteCallback` is used to detect the jump times.
When `jump_interp_pts` is set to a positive integer, a `ContinuousCallback` is used to detect the jump times.
"""
function mcsolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector, c_ops::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    jump_interp_pts::Integer=-1,
    kwargs...) where {T1,T2}

    H_eff = H - T2(0.5im) * mapreduce(op -> op' * op, +, c_ops)

    cb1 = jump_interp_pts == -1 ? DiscreteLindbladJumpCallback() : ContinuousLindbladJumpCallback(jump_interp_pts)
    kwargs2 = kwargs
    kwargs2 = merge(kwargs2, haskey(kwargs2, :callback) ? 
                Dict(:callback => CallbackSet(cb1, kwargs2[:callback])) : Dict(:callback => cb1))
    if !isempty(e_ops)
        cb2 = PresetTimeCallback(t_l, _save_func_mcsolve, save_positions=(false, false))
        kwargs2 = merge(kwargs2, Dict(:callback => CallbackSet(kwargs2[:callback], cb2)))
    end

    e_ops2 = length(e_ops) == 0 ? Vector{Matrix{T1}}([]) : get_data.(e_ops)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    cache_mc = similar(ψ0.data)
    weights_mc = Array{Float64}(undef, length(c_ops))
    cumsum_weights_mc = similar(weights_mc)
    params2 = merge(params, Dict(:expvals => expvals, :e_ops => e_ops2, :save_it => Ref{Int32}(0), 
                                :random_n => Ref(rand()), :c_ops => get_data.(c_ops), :cache_mc => cache_mc, 
                                :weights_mc => weights_mc, :cumsum_weights_mc => cumsum_weights_mc,
                                :jump_times => Float64[], :jump_which => Int16[]))

    return sesolveProblem(H_eff, ψ0, t_l; alg=alg, H_t=H_t, params=params2, progress=false, kwargs2...)
end

"""
    mcsolveEnsembleProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        n_traj::Integer=1,
        jump_interp_pts::Integer=-1,
        prob_func::Function=_mcsolve_prob_func,
        output_func::Function=_mcsolve_output_func,
        kwargs...)

Generates the ODEProblem for an ensemble of trajectories of the Monte Carlo wave function
time evolution of an open quantum system.

# Arguments
- `H::QuantumObject`: Hamiltonian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::AbstractVector`: List of collapse operators.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function}`: Time-dependent part of the Hamiltonian.
- `params::Dict{Symbol, Any}`: Dictionary of parameters to pass to the solver.
- `progress::Bool`: Whether to show a progress bar.
- `n_traj::Integer`: Number of trajectories to use.
- `jump_interp_pts::Integer`: Number of points to use for interpolation of the jump times.
- `prob_func::Function`: Function to use for generating the ODEProblem.
- `output_func::Function`: Function to use for generating the output of a single trajectory.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns
- `prob::EnsembleProblem with ODEProblem`: The Ensemble ODEProblem for the Monte Carlo
wave function time evolution.

# Notes
When `jump_interp_pts` is set to -1, a `DiscreteCallback` is used to detect the jump times.
When `jump_interp_pts` is set to a positive integer, a `ContinuousCallback` is used to detect the jump times.
"""
function mcsolveEnsembleProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector, c_ops::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    n_traj::Integer=1,
    jump_interp_pts::Integer=-1,
    prob_func::Function=_mcsolve_prob_func,
    output_func::Function=_mcsolve_output_func,
    kwargs...) where {T1,T2}

    progr = Progress(n_traj, showspeed=true, enabled=progress)
    channel = RemoteChannel(() -> Channel{Bool}(), 1)
    @async while take!(channel)
        next!(progr)
    end

    params2 = merge(params, Dict(:progr_channel => channel))

    prob_mc = mcsolveProblem(H, ψ0, t_l, c_ops; alg=alg, e_ops=e_ops, 
                H_t=H_t, params=params2, jump_interp_pts=jump_interp_pts, kwargs...)


    ensemble_prob = EnsembleProblem(prob_mc, prob_func=prob_func,
                            output_func=output_func, safetycopy=false)

    return ensemble_prob
end


"""
    sesolve(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        kwargs...)

Time evolution of a closed quantum system using the Schrödinger equation.

# Arguments
- `H::QuantumObject`: Hamiltonian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function}`: Time-dependent part of the Hamiltonian.
- `params::Dict{Symbol, Any}`: Dictionary of parameters to pass to the solver.
- `progress::Bool`: Whether to show a progress bar.
- `kwargs...`: Additional keyword arguments to pass to the solver.

- Returns
- `sol::TimeEvolutionSol`: The solution of the time evolution.
"""
function sesolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    kwargs...) where {T1,T2}

    prob = sesolveProblem(H, ψ0, t_l; alg=alg, e_ops=e_ops,
            H_t=H_t, params=params, progress=progress, kwargs...)
    
    return sesolve(prob; alg=alg, kwargs...)
end

function sesolve(prob::ODEProblem;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    kwargs...)

    sol = solve(prob, alg)

    Hdims = sol.prob.p.Hdims
    ψt = !haskey(kwargs, :save_idxs) ? map(ϕ -> QuantumObject(ϕ, dims = Hdims), sol.u) : sol.u

    return TimeEvolutionSol(sol.t, ψt, sol.prob.p.expvals)
end

"""
    mesolve(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector=[];
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        kwargs...)

Time evolution of an open quantum system using master equation.

# Arguments
- `H::QuantumObject`: Hamiltonian of Liouvillian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::AbstractVector`: List of collapse operators.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function}`: Time-dependent part of the Hamiltonian.
- `params::Dict{Symbol, Any}`: Dictionary of parameters to pass to the solver.
- `progress::Bool`: Whether to show a progress bar.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns
- `sol::TimeEvolutionSol`: The solution of the time evolution.
"""
function mesolve(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector=[];
    alg::OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    kwargs...) where {T1,T2,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    prob = mesolveProblem(H, ψ0, t_l, c_ops; alg=alg, e_ops=e_ops,
            H_t=H_t, params=params, progress=progress, kwargs...)
    
    return mesolve(prob; alg=alg, kwargs...)
end

function mesolve(prob::ODEProblem;
    alg::OrdinaryDiffEqAlgorithm=Tsit5(),
    kwargs...)

    sol = solve(prob, alg)

    Hdims = sol.prob.p.Hdims
    ρt = !haskey(kwargs, :save_idxs) ? map(ϕ -> QuantumObject(vec2mat(ϕ), dims=Hdims), sol.u) : sol.u

    return TimeEvolutionSol(sol.t, ρt, sol.prob.p.expvals)
end

"""
    mcsolve(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        n_traj::Integer=1,
        ensemble_method=EnsembleThreads(),
        jump_interp_pts::Integer=-1,
        kwargs...)

Time evolution of an open quantum system using quantum trajectories.

# Arguments
- `H::QuantumObject`: Hamiltonian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `c_ops::AbstractVector`: List of collapse operators.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function}`: Time-dependent part of the Hamiltonian.
- `params::Dict{Symbol, Any}`: Dictionary of parameters to pass to the solver.
- `progress::Bool`: Whether to show a progress bar.
- `n_traj::Integer`: Number of trajectories to use.
- `ensemble_method`: Ensemble method to use.
- `jump_interp_pts::Integer`: Number of points to use for interpolation of jump times.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns
- `sol::TimeEvolutionMCSol`: The solution of the time evolution.

# Notes
`ensemble_method` can be one of `EnsembleThreads()`, `EnsembleSerial()`, `EnsembleDistributed()`.
When `jump_interp_pts` is set to -1, a `DiscreteCallback` is used to detect the jump times.
When `jump_interp_pts` is set to a positive integer, a `ContinuousCallback` is used to detect the jump times.
"""
function mcsolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector, c_ops::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    n_traj::Integer=1,
    ensemble_method=EnsembleThreads(),
    jump_interp_pts::Integer=-1,
    kwargs...) where {T1,T2}


    ens_prob_mc = mcsolveEnsembleProblem(H, ψ0, t_l, c_ops; alg=alg, e_ops=e_ops, 
                H_t=H_t, params=params, progress=progress, n_traj=n_traj, 
                jump_interp_pts=jump_interp_pts, kwargs...)

    return mcsolve(ens_prob_mc; alg=alg, n_traj=n_traj, ensemble_method=ensemble_method, kwargs...)
end

function mcsolve(ens_prob_mc::EnsembleProblem;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
    n_traj::Integer=1,
    ensemble_method=EnsembleThreads(),
    kwargs...)

    sol = solve(ens_prob_mc, alg, ensemble_method, trajectories=n_traj)
    put!(sol[1].prob.p.progr_channel, false)

    expvals_all = Array{ComplexF64}(undef, length(sol), size(sol[1].prob.p.expvals)...)
    times = Vector{Vector{Float64}}([])
    states = haskey(sol[1].prob.kwargs, :save_idxs) ? Vector{Vector{eltype(sol[1].u[1])}}([]) : Vector{Vector{QuantumObject}}([])
    jump_times = Vector{Vector{Float64}}([])
    jump_which = Vector{Vector{Int16}}([])
    foreach(i -> _mcsolve_generate_statistics(sol, i, times, states, expvals_all, jump_times, jump_which), eachindex(sol))
    expvals = dropdims(sum(expvals_all, dims=1), dims=1) ./ length(sol)

    TimeEvolutionMCSol(times, states, expvals, expvals_all, jump_times, jump_which)
end








### LIOUVILLIAN AND STEADYSTATE ###
@doc raw"""
    liouvillian(H::QuantumObject, c_ops::AbstractVector)

Construct the Liouvillian superoperator for a system Hamiltonian and a set of collapse operators:
``\mathcal{L} \cdot = -i[\hat{H}, \cdot] + \sum_i \mathcal{D}[\hat{O}_i] \cdot``,
where ``\mathcal{D}[\hat{O}_i] \cdot = \hat{O}_i \cdot \hat{O}_i^\dagger - \frac{1}{2} \hat{O}_i^\dagger \hat{O}_i \cdot - \frac{1}{2} \cdot \hat{O}_i^\dagger \hat{O}_i``.
"""
function liouvillian(H::QuantumObject{<:AbstractArray{T},OpType},
    c_ops::AbstractVector) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}

    L = isoper(H) ? -1im * (spre(H) - spost(H)) : H
    for c_op in c_ops
        isoper(c_op) ? L += lindblad_dissipator(c_op) : L += c_op
    end
    L
end


liouvillian(H::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    liouvillian(H, [])

function liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real; n_max::Int=4, solver::LSolver=LiouvillianDirectSolver()) where {T1,T2,T3,LSolver<:LiouvillianSolver}

    ((L₀.dims == Lₚ.dims) && (L₀.dims == Lₘ.dims)) || throw(ErrorException("The operators are not of the same Hilbert dimension."))

    _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, solver, n_max=n_max)
end

function liouvillian_floquet(H::QuantumObject{<:AbstractArray{T1},OpType1},
    c_ops::AbstractVector,
    Hₚ::QuantumObject{<:AbstractArray{T2},OpType2},
    Hₘ::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, solver::LSolver=LiouvillianDirectSolver()) where {T1,T2,T3,
                                                                            OpType1<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            OpType2<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            OpType3<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                            LSolver<:LiouvillianSolver}

    liouvillian_floquet(liouvillian(H, c_ops), liouvillian(Hₚ), liouvillian(Hₘ), ω, solver=solver, n_max=n_max)
end

@doc raw"""
    liouvillian_generalized(H::QuantumObject, fields::Vector, 
    κ_list::Vector, ω_list::Vector, T_list::Vector; N_trunc::Int=size(H,1), tol::Float64=0.0)

Constructs the generalized Liouvillian for a system coupled to a bath of harmonic oscillators.

See, e.g., Settineri, Alessio, et al. "Dissipation and thermal noise in hybrid quantum systems in the ultrastrong-coupling regime." Physical Review A 98.5 (2018): 053834.
"""
function QuPhys.liouvillian_generalized(H::QuantumObject{<:AbstractArray, OperatorQuantumObject}, fields::Vector, 
    κ_list::Vector{<:Number}, ω_list::Vector{<:Number}, T_list::Vector{<:Real}; N_trunc::Int=size(H,1), tol::Real=1e-14)

    (length(fields) == length(κ_list) == length(ω_list) == length(T_list)) || throw(DimensionMismatch("The number of fields, κs, ωs and Ts must be the same."))

    dims = N_trunc == size(H,1) ? H.dims : [N_trunc]
    E, U = eigen(H)
    E = E[1:N_trunc]
    U = QuantumObject(U, dims=H.dims)

    H_d = QuantumObject(spdiagm(complex(E)), dims=dims)

    Ω = dense_to_sparse(E' .- E, tol)
    Ω = triu(Ω, 1)

    L = liouvillian(H_d)

    for i in eachindex(fields)
        X_op = dense_to_sparse((U' * fields[i] * U).data[1:N_trunc, 1:N_trunc], tol)
        if ishermitian(fields[i])
            X_op = (X_op + X_op') / 2 # Make sure it's hermitian
        end

        # Nikki Master Equation
        N_th = n_th.(Ω, T_list[i])
        Sp₀ = QuantumObject( triu(X_op, 1), dims=dims )
        Sp₁ = QuantumObject( droptol!( (@. (Ω / ω_list[i]) * N_th * Sp₀.data), tol), dims=dims )
        Sp₂ = QuantumObject( droptol!( (@. (Ω / ω_list[i]) * (1 + N_th) * Sp₀.data), tol), dims=dims )
        S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

        L += κ_list[i] / 2 * ( sprepost(Sp₁', Sp₀) + sprepost(Sp₀', Sp₁) - spre(Sp₀ * Sp₁') - spost(Sp₁ * Sp₀') )
        L += κ_list[i] / 2 * ( sprepost(Sp₂, Sp₀') + sprepost(Sp₀, Sp₂') - spre(Sp₀' * Sp₂) - spost(Sp₂' * Sp₀) )
        L += κ_list[i] * T_list[i] / (4 * ω_list[i]) * ( 4 * sprepost(S0, S0) - 2 * spre(S0 * S0) - 2 * spost(S0 * S0) )
    end

    return E, U, L
end

function _liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real, solver::LiouvillianDirectSolver; n_max::Int=4) where {T1,T2,T3}

    L_0 = L₀.data
    L_p = sparse_to_dense(Lₚ.data)
    L_m = sparse_to_dense(Lₘ.data)

    n_i = n_max
    S = -(L_0 - 1im * n_i * ω * I) \ L_p
    T = -(L_0 + 1im * n_i * ω * I) \ L_m

    for n_i in n_max-1:-1:1
        S = -(L_0 - 1im * n_i * ω * I + L_m * S) \ L_p
        T = -(L_0 + 1im * n_i * ω * I + L_p * T) \ L_m
    end

    solver.tol == 0 && return QuantumObject(L_0 + L_m * S + L_p * T, SuperOperatorQuantumObject, L₀.dims)
    return QuantumObject(dense_to_sparse(L_0 + L_m * S + L_p * T, solver.tol), SuperOperatorQuantumObject, L₀.dims)
end

function steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject};
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    _steadystate(L, solver)
end

function steadystate(H::QuantumObject{<:AbstractArray{T},OpType}, c_ops::Vector,
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},SSSolver<:SteadyStateSolver}

    L = liouvillian(H, c_ops)
    steadystate(L, solver=solver)
end

function _steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::Type{SteadyStateDirectSolver}) where {T}

    L_tmp = copy(L.data)
    N = prod(L.dims)
    weight = norm(L_tmp, 1) / length(L_tmp)
    v0 = zeros(ComplexF64, N^2) # This is not scalable for GPU arrays
    v0[1] = weight

    L_tmp[1, [N * (i - 1) + i for i in 1:N]] .+= weight

    ρss_vec = L_tmp \ v0
    ρss = reshape(ρss_vec, N, N)
    ρss = (ρss + ρss') / 2 # Hermitianize
    QuantumObject(ρss, OperatorQuantumObject, L.dims)
end

@doc raw"""
    steadystate_floquet(H_0::QuantumObject,
        c_ops::Vector, H_p::QuantumObject,
        H_m::QuantumObject,
        ω::Real; n_max::Int=4, lf_solver::LSolver=LiouvillianDirectSolver(),
        ss_solver::Type{SSSolver}=SteadyStateDirectSolver)

Calculates the steady state of a periodically driven system.
Here `H_0` is the Hamiltonian or the Liouvillian of the undriven system.
Considering a monochromatic drive at frequency ``\\omega``, we divide it into two parts,
`H_p` and `H_m`, where `H_p` oscillates
as ``e^{i \\omega t}`` and `H_m` oscillates as ``e^{-i \\omega t}``.
`n_max` is the number of iterations used to obtain the effective Liouvillian,
`lf_solver` is the solver used to solve the effective Liouvillian,
and `ss_solver` is the solver used to solve the steady state.
"""
function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T1},OpType1},
    c_ops::Vector, H_p::QuantumObject{<:AbstractArray{T2},OpType2},
    H_m::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, lf_solver::LSolver=LiouvillianDirectSolver(),
    ss_solver::Type{SSSolver}=SteadyStateDirectSolver) where {T1,T2,T3,OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    LSolver<:LiouvillianSolver,SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end
