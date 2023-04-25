abstract type LiouvillianSolver end
abstract type LiouvillianDirectSolver <: LiouvillianSolver end

abstract type SteadyStateSolver end
abstract type SteadyStateDirectSolver <: SteadyStateSolver end

struct TimeEvolutionSol
    times::AbstractVector
    states::AbstractArray
    expect::AbstractArray
end

function _save_func_sesolve(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params["progr"]
    e_ops = internal_params["e_ops"]
    expvals = internal_params["expvals"]
    expvals[:, progr.counter+1] .= map(op -> dot(u, op.data, u), e_ops)
    next!(progr)
end

function _save_func_mesolve(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params["progr"]
    e_ops = internal_params["e_ops"]
    expvals = internal_params["expvals"]
    expvals[:, progr.counter+1] .= map(op -> tr(op.data * reshape(u, isqrt(length(u)), :)), e_ops)
    next!(progr)
end

function _save_func_mcsolve(u, t, integrator)
    internal_params = integrator.p[1]
    save_it = internal_params["save_it"]
    e_ops = internal_params["e_ops"]
    expvals = internal_params["expvals"]
    ψ = normalize(u)
    expvals[:, save_it[]+1] .= map(op -> dot(ψ, op.data, ψ), e_ops)
    save_it[]+=1
end

function LindbladJumpAffect!(integrator)
    ψ = integrator.u
    internal_params = integrator.p[1]
    c_ops = internal_params["c_ops"]

    if length(c_ops) == 1
        integrator.u = normalize!(c_ops[1] * ψ)
    else
        collaps_idx = 1
        r2 = rand()
        dp = 0
        c_op_ψ_l = Vector{Float64}(undef, length(c_ops))
        @inbounds for i in eachindex(c_ops)
            c_op_ψ = c_ops[i] * ψ
            res = real(dot(c_op_ψ, c_op_ψ))
            c_op_ψ_l[i] = res
            dp += res
        end
        prob = 0
        @inbounds for i in eachindex(c_ops)
            res = c_op_ψ_l[i]
            prob += res / dp
            if prob >= r2
                collaps_idx = i
                break
            end
        end
        integrator.u = normalize!(c_ops[collaps_idx] * ψ)
    end
    integrator.p[1]["random_n"] = rand()
end

function ContinuousLindbladJumpCallback(interp_points::Int=0; affect_func::Function=LindbladJumpAffect!)
    LindbladJumpCondition(u, t, integrator) = integrator.p[1]["random_n"] - norm(u)^2

    ContinuousCallback(LindbladJumpCondition, affect_func, nothing, interp_points=interp_points, save_positions=(false, false))
end

function DiscreteLindbladJumpCallback(;affect_func::Function=LindbladJumpAffect!)
    LindbladJumpCondition(u, t, integrator) = norm(u)^2 < integrator.p[1]["random_n"]

    DiscreteCallback(LindbladJumpCondition, affect_func, save_positions=(false, false))
end

"""
    sesolve(H::QuantumObject,
                ψ0::QuantumObject,
                t_l::AbstractVector;  
                e_ops::AbstractVector=[], 
                alg=Vern7(),
                H_t=nothing, 
                params::AbstractVector=[],
                progress::Bool=true,
                callbacks=[],
                kwargs...)

Time evolution of a closed quantum system using the Schrödinger equation.
"""
function sesolve(H::QuantumObject{<:AbstractArray{T},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
    t_l::AbstractVector;
    e_ops::AbstractVector=[],
    alg=Vern7(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks=[],
    kwargs...) where {T}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims

    tspan = (t_l[1], t_l[end])

    H0 = -1im * H.data
    ψ0 = ψ0.data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    p = [Dict("H" => H0, "e_ops" => e_ops, "expvals" => expvals, "progr" => progr), params...]

    is_time_dependent = !(H_t === nothing)

    cb1 = FunctionCallingCallback(_save_func_sesolve, funcat=t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, callbacks...)

    if typeof(alg) <: LinearExponential
        is_time_dependent && error("The Hamiltonian must to be time independent when using LinearExponential algorithm.")
        A = DiffEqArrayOperator(H0)
        prob = ODEProblem(A, ψ0, tspan, p; kwargs...)
        sol = solve(prob, alg, dt=(t_l[2] - t_l[1]), callback=cb)
    else
        if !is_time_dependent
            dudt! = (du, u, p, t) -> mul!(du, p[1]["H"], u)
        else
            dudt! = (du, u, p, t) -> mul!(du, p[1]["H"] - 1im * H_t(t).data, u)
        end
        prob = ODEProblem(dudt!, ψ0, tspan, p; kwargs...)
        sol = solve(prob, alg, callback=cb)
    end

    ψt_len = length(sol.u[1])
    ψt_len == prod(Hdims) ? ψt = [QuantumObject(normalize!(ϕ), dims=Hdims) for ϕ in sol.u] : ψt = sol.u

    return TimeEvolutionSol(sol.t, ψt, sol.prob.p[1]["expvals"])
end

"""
    mesolve(H::QuantumObject,
            ψ0::QuantumObject,
            t_l::AbstractVector, c_ops::AbstractVector=[]; 
            e_ops::AbstractVector=[], 
            alg=Vern7(), 
            H_t=nothing, 
            params::AbstractVector=[],
            progress::Bool = true,
            callbacks = [],
            kwargs...)

Time evolution of an open quantum system using master equation.
"""
function mesolve(H::QuantumObject{<:AbstractArray{T},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector=[];
    e_ops::AbstractVector=[],
    alg=Vern7(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks=[],
    kwargs...) where {T,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ρ0 = isket(ψ0) ? reshape(ψ0.data * ψ0.data', length(ψ0)^2) : reshape(ψ0.data, length(ψ0))

    L = liouvillian(H, c_ops).data

    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    p = [Dict("L" => L, "e_ops" => e_ops, "expvals" => expvals, "progr" => progr), params...]

    is_time_dependent = !(H_t === nothing)

    cb1 = FunctionCallingCallback(_save_func_mesolve, funcat=t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, callbacks...)

    if typeof(alg) <: LinearExponential
        is_time_dependent && error("The Liouvillian must to be time independent when using LinearExponential algorith.")
        A = DiffEqArrayOperator(L)
        prob = ODEProblem(A, ρ0, tspan, p; kwargs...)
        sol = solve(prob, alg, dt=(t_l[2] - t_l[1]), callback=cb)
    else
        if !is_time_dependent
            dudt! = (du, u, p, t) -> mul!(du, p[1]["L"], u)
        else
            if H_t(0.0).type <: OperatorQuantumObject
                @warn string("To speed up the calculation, it is always better to define ",
                    "the time-dependent part as a SuperOperator, and not as an Operator.") maxlog = 1
                dudt! = (du, u, p, t) -> mul!(du, p[1]["L"] + liouvillian(H_t(t)).data, u)
            else
                dudt! = (du, u, p, t) -> mul!(du, p[1]["L"] + H_t(t).data, u)
            end
        end
        prob = ODEProblem(dudt!, ρ0, tspan, p; kwargs...)
        sol = solve(prob, alg, callback=cb)
    end

    ρt_len = isqrt(length(sol.u[1]))
    ρt_len == prod(Hdims) ? ρt = [QuantumObject(reshape(ϕ, ρt_len, ρt_len), dims=Hdims) for ϕ in sol.u] : ρt = sol.u

    return TimeEvolutionSol(sol.t, ρt, sol.prob.p[1]["expvals"])
end

"""
    mcsolve(H::QuantumObject,
            ψ0::QuantumObject,
            t_l::AbstractVector, c_ops::AbstractVector;
            e_ops::AbstractVector=[],
            n_traj::Int=1,
            batch_size::Int=min(Threads.nthreads(), n_traj),
            alg=AutoVern7(KenCarp4(autodiff=false)),
            ensemble_method=EnsembleThreads(),
            H_t=nothing,
            params::AbstractVector=[],
            progress::Bool=true,
            jump_interp_pts::Int=10,
            callbacks=[],
            kwargs...)

Time evolution of an open quantum system using quantum trajectories.
"""
function mcsolve(H::QuantumObject{<:AbstractArray{T},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
    t_l::AbstractVector, c_ops::AbstractVector;
    e_ops::AbstractVector=[],
    n_traj::Int=1,
    batch_size::Int=min(Threads.nthreads(), n_traj),
    alg=AutoVern7(KenCarp4(autodiff=false)),
    ensemble_method=EnsembleThreads(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    jump_interp_pts::Int=10,
    callbacks=[],
    kwargs...) where {T}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims

    tspan = (t_l[1], t_l[end])
    e_ops_len = length(e_ops)

    ψ0 = ψ0.data

    progr = Progress(n_traj, showspeed=true, enabled=progress)
    channel = RemoteChannel(() -> Channel{Bool}(), 1)
    @async while take!(channel)
        next!(progr)
    end

    function prob_func(prob, i, repeat)
        H_eff = copy(H.data)
        for c_op in c_ops
            H_eff += -0.5im * c_op.data' * c_op.data
        end
        H_eff = -1im * H_eff
        c_ops0 = map(op -> op.data, c_ops)
        
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        p = [Dict("H" => H_eff, "c_ops" => c_ops0, "e_ops" => e_ops, "random_n" => rand(), "expvals" => expvals, "save_it" => Ref{Int32}(0)), params...]
        remake(prob, p=p)
    end
    function output_func(sol, i)
        put!(channel, true)
        if e_ops_len == 0
            res = map(ϕ -> QuantumObject(ϕ, dims=Hdims), sol.u)
        else
            res = sol.prob.p[1]["expvals"]
        end
        (res, false)
    end
    function reduction(u, batch, I)
        if e_ops_len == 0
            tmp = hcat(batch...)
            length(u) == 0 && return tmp, false
            res = hcat(u, tmp)
        else
            tmp = sum(cat(batch..., dims=3), dims=3)
            length(u) == 0 && return tmp, false
            res = sum(cat(u, tmp, dims=3), dims=3)
        end
        return res, false
    end

    is_time_dependent = !(H_t === nothing)
    if is_time_dependent
        dudt! = (du, u, p, t) -> mul!(du, p[1]["H"] - 1im * H_t(t).data, u)
    else
        dudt! = (du, u, p, t) -> mul!(du, p[1]["H"], u)
    end

    cb1 = FunctionCallingCallback(_save_func_mcsolve, funcat=t_l)
    cb2 = jump_interp_pts == -1 ? DiscreteLindbladJumpCallback() : ContinuousLindbladJumpCallback(jump_interp_pts)
    cb3 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, cb3, callbacks...)

    prob = ODEProblem(dudt!, ψ0, tspan, callback=cb; kwargs...)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, reduction=reduction)
    sol = solve(ensemble_prob, alg, ensemble_method, trajectories=n_traj, batch_size=batch_size)

    put!(channel, false)

    e_ops_len == 0 && return TimeEvolutionSol(t_l, sol.u, [])

    e_ops_expect = sum(sol.u, dims=3) ./ n_traj

    return TimeEvolutionSol(t_l, [], e_ops_expect)
end


### LIOUVILLIAN AND STEADYSTATE ###
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
    Lₚ::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    ω::Real; n_max::Int=4, solver::Type{LSolver}=LiouvillianDirectSolver) where {T1,LSolver<:LiouvillianSolver}

    ((L₀.dims == Lₚ.dims) && (L₀.dims == Lₘ.dims)) || throw(ErrorException("The operators are not of the same Hilbert dimension."))

    _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, solver, n_max=n_max)
end

function _liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    ω::Real, solver::Type{LiouvillianDirectSolver}; n_max::Int=4) where {T1}

    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data
    S = T = spzeros(T1, size(L_0)...)

    L_p_d = Matrix(L_p)
    L_m_d = Matrix(L_m)

    for n_i in n_max:-1:1
        S, T = -(L_0 - 1im * n_i * ω * I + L_m_d * S) \ L_p_d, -(L_0 + 1im * n_i * ω * I + L_p_d * T) \ L_m_d
    end

    QuantumObject(L_0 + L_m * S + L_p * T, SuperOperatorQuantumObject, L₀.dims)
end

function steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject};
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    _steadystate(L, solver)
end

function steadystate(H::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, c_ops::Vector,
    solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    L = liouvillian(H, c_ops)
    steadystate(L, solver=solver)
end

function _steadystate(L::QuantumObject{<:AbstractArray{T},SuperOperatorQuantumObject},
    solver::Type{SteadyStateDirectSolver}) where {T}

    L_tmp = copy(L.data)
    N = prod(L.dims) # floor(Int, √(size(L_tmp, 1)))
    weight = sum(abs.(L_tmp)) / length(L_tmp)
    v0 = zeros(ComplexF64, N^2)
    v0[1] = weight

    L_tmp[1, [N * (i - 1) + i for i in 1:N]] .+= weight

    rho_ss_vec = L_tmp \ v0
    rho_ss = reshape(rho_ss_vec, N, N)
    QuantumObject(rho_ss, OperatorQuantumObject, L.dims)
end

function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T},OpType1},
    c_ops::Vector, H_p::QuantumObject{<:AbstractArray{T},OpType2},
    H_m::QuantumObject{<:AbstractArray{T},OpType3},
    ω::Real; n_max::Int=4, lf_solver::Type{LSolver}=LiouvillianDirectSolver,
    ss_solver::Type{SSSolver}=SteadyStateDirectSolver) where {T,OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    LSolver<:LiouvillianSolver,SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end
