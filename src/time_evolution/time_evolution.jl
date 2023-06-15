abstract type LiouvillianSolver end
struct LiouvillianDirectSolver <: LiouvillianSolver 
    tol::Real
end

abstract type SteadyStateSolver end
abstract type SteadyStateDirectSolver <: SteadyStateSolver end

struct TimeEvolutionSol
    times::AbstractVector
    states::AbstractArray
    expect::AbstractArray
end

LiouvillianDirectSolver(;tol=1e-14) = LiouvillianDirectSolver(tol)

mutable struct SesolveParams{LT<:AbstractArray, TE<:AbstractArray, T<:BlasFloat}
    U::LT # Liouvillian
    e_ops::Vector{TE} # expectation operators
    expvals::Matrix{T} # expectation values
    progr::Progress # progress bar
end

mutable struct MesolveParams{LT<:AbstractArray, TE<:AbstractArray, T<:BlasFloat}
    L::LT # Liouvillian
    e_ops::Vector{TE} # expectation operators
    expvals::Matrix{T} # expectation values
    progr::Progress # progress bar
end

function _save_func_sesolve(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params.progr
    e_ops = internal_params.e_ops
    expvals = internal_params.expvals

    if !isempty(e_ops)
        _expect = op -> dot(u, op, u)
        @. expvals[:, progr.counter+1] = _expect(e_ops)
    end
    next!(progr)
end

function _save_func_mesolve(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params.progr
    e_ops = internal_params.e_ops
    expvals = internal_params.expvals

    if !isempty(e_ops)
        # This is equivalent to tr(op * ρ), when both are matrices.
        # The advantage of using this convention is that I don't need
        # to reshape u to make it a matrix, but I reshape the e_ops once.
        # Attention! dot(u, op) is different from dot(op, u), since this function
        # automatically conjugates the left vector, which leads to wrong results for
        # non-hermitian operators.
        
        _expect = op -> dot(op, u)
        @. expvals[:, progr.counter+1] = _expect(e_ops)
    end
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
    sesolveProblem(H::QuantumObject,
    ψ0::QuantumObject,
    t_l::AbstractVector,
    alg;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...)

Generates the ODEProblem for the Schrödinger time evolution of a quantum system.
"""
function sesolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ϕ0 = get_data(ψ0)
    U = -1im * H.data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = get_data.(e_ops)
    p = [SesolveParams(U, e_ops2, expvals, progr), params...]

    cb1 = FunctionCallingCallback(_save_func_sesolve, funcat=t_l)
    if haskey(kwargs, :callback)
        kwargs2 = merge(kwargs, Dict(:callback => CallbackSet(cb1, kwargs[:callback])))
    else
        kwargs2 = merge(kwargs, Dict(:callback => cb1))
    end
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])
    if !is_time_dependent
        dudt! = (du, u, p, t) -> mul!(du, p[1].U, u)
    else
        dudt! = (du, u, p, t) -> mul!(du, p[1].U - 1im * H_t(t).data, u)
    end

    ODEProblem(dudt!, ϕ0, tspan, p; kwargs2...)
end

function sesolveProblem(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqExponentialAlgorithm;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)
    is_time_dependent && error("The Liouvillian must to be time independent when using LinearExponential algorith.")

    ϕ0 = get_data(ψ0)
    U = -1im * H.data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = get_data.(e_ops)
    p = [SesolveParams(U, e_ops2, expvals, progr), params...]

    cb1 = FunctionCallingCallback(_save_func_sesolve, funcat=t_l)
    if haskey(kwargs, :callback)
        kwargs2 = merge(kwargs, Dict(:callback => CallbackSet(cb1, kwargs[:callback])))
    else
        kwargs2 = merge(kwargs, Dict(:callback => cb1))
    end
    kwargs2 = merge(kwargs2, Dict(:dt => t_l[2] - t_l[1]))
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])
    A = DiffEqArrayOperator(U)
    ODEProblem(A, ϕ0, tspan, p; kwargs2...)
end

"""
    mesolveProblem(H::QuantumObject,
    ψ0::QuantumObject,
    t_l::AbstractVector, c_ops::AbstractVector,
    alg;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...)

Generates the ODEProblem for the master equation time evolution of an open quantum system.
"""
function mesolveProblem(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ρ0 = isket(ψ0) ? mat2vec(ket2dm(ψ0).data) : mat2vec(ψ0.data)
    L = liouvillian(H, c_ops).data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = isempty(e_ops) ? Vector{Vector}([]) : map(op -> mat2vec(sparse_to_dense(get_data(op)')), e_ops)
    p = [MesolveParams(L, e_ops2, expvals, progr), params...]

    cb1 = FunctionCallingCallback(_save_func_mesolve, funcat=t_l)
    if haskey(kwargs, :callback)
        kwargs2 = merge(kwargs, Dict(:callback => CallbackSet(cb1, kwargs[:callback])))
    else
        kwargs2 = merge(kwargs, Dict(:callback => cb1))
    end
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])
    if !is_time_dependent
        dudt! = (du, u, p, t) -> mul!(du, p[1].L, u)
    else
        if isoper(H_t(0.0))
            @warn string("To speed up the calculation, it is always better to define ",
                "the time-dependent part as a SuperOperator, and not as an Operator.") maxlog = 1
            dudt! = (du, u, p, t) -> mul!(du, p[1].L + liouvillian(H_t(t)).data, u)
        else
            dudt! = (du, u, p, t) -> mul!(du, p[1].L + H_t(t).data, u)
        end
    end
    
    ODEProblem(dudt!, ρ0, tspan, p; kwargs2...)
end

function mesolveProblem(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector,
    alg::OrdinaryDiffEq.OrdinaryDiffEqExponentialAlgorithm;
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)
    is_time_dependent && error("The Liouvillian must to be time independent when using LinearExponential algorith.")

    ρ0 = isket(ψ0) ? mat2vec(ket2dm(ψ0).data) : mat2vec(ψ0.data)
    L = liouvillian(H, c_ops).data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = map(op -> mat2vec(sparse_to_dense(get_data(op))), e_ops)
    p = [MesolveParams(L, e_ops2, expvals, progr), params...]

    cb1 = FunctionCallingCallback(_save_func_mesolve, funcat=t_l)
    if haskey(kwargs, :callback)
        kwargs2 = merge(kwargs, Dict(:callback => CallbackSet(cb1, kwargs[:callback])))
    else
        kwargs2 = merge(kwargs, Dict(:callback => cb1))
    end
    kwargs2 = merge(kwargs2, Dict(:dt => t_l[2] - t_l[1]))
    !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, Dict(:abstol => 1e-7)))
    !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, Dict(:reltol => 1e-5)))
    !haskey(kwargs2, :saveat) && (kwargs2 = merge(kwargs2, Dict(:saveat => [t_l[end]])))

    tspan = (t_l[1], t_l[end])
    A = DiffEqArrayOperator(L)
    ODEProblem(A, ρ0, tspan, p; kwargs2...)
end




"""
    sesolve(H::QuantumObject,
                ψ0::QuantumObject,
                t_l::AbstractVector;
                alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
                e_ops::AbstractVector=[], 
                H_t=nothing, 
                params::AbstractVector=[],
                progress::Bool=true,
                callbacks=[],
                kwargs...)

Time evolution of a closed quantum system using the Schrödinger equation.
"""
function sesolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2}

    prob = sesolveProblem(H, ψ0, t_l, alg; e_ops=e_ops,
            H_t=H_t, params=params, progress=progress, kwargs...)
    sol = solve(prob, alg)

    ψt = length(sol.u[1]) == prod(H.dims) ? map(ϕ -> QuantumObject(ϕ, dims=H.dims), sol.u) : sol.u

    return TimeEvolutionSol(sol.t, ψt, sol.prob.p[1].expvals)
end

"""
    mesolve(H::QuantumObject,
            ψ0::QuantumObject,
            t_l::AbstractVector, c_ops::AbstractVector=[]; 
            alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
            e_ops::AbstractVector=[], 
            H_t=nothing, 
            params::AbstractVector=[],
            progress::Bool = true,
            callbacks = [],
            kwargs...)

Time evolution of an open quantum system using master equation.
"""
function mesolve(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, c_ops::AbstractVector=[];
    alg::OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::AbstractVector=[],
    H_t::Union{Nothing,Function}=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    kwargs...) where {T1,T2,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    prob = mesolveProblem(H, ψ0, t_l, c_ops, alg; e_ops=e_ops,
            H_t=H_t, params=params, progress=progress, kwargs...)
    sol = solve(prob, alg)

    ρt = length(sol.u[1]) == prod(H.dims)^2 ? map(ϕ -> QuantumObject(reshape(ϕ, prod(H.dims), prod(H.dims)), dims=H.dims), sol.u) : sol.u

    return TimeEvolutionSol(sol.t, ρt, sol.prob.p[1].expvals)
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
function mcsolve(H::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
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
    kwargs...) where {T1,T2}

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

@doc raw"""
    liouvillian_generalized(H::QuantumObject, fields::Vector, 
    κ_list::Vector, ω_list::Vector, T_list::Vector; N_trunc::Int=size(H,1), tol::Float64=0.0)

Constructs the generalized Liouvillian for a system coupled to a bath of harmonic oscillators.

See, e.g., Settineri, Alessio, et al. "Dissipation and thermal noise in hybrid quantum systems in the ultrastrong-coupling regime." Physical Review A 98.5 (2018): 053834.
"""
function liouvillian_generalized(H::QuantumObject{<:AbstractArray, OperatorQuantumObject}, fields::Vector, 
    κ_list::Vector{<:Number}, ω_list::Vector{<:Number}, T_list::Vector{<:Real}; N_trunc::Int=size(H,1), tol::Real=1e-14)

    (length(fields) == length(κ_list) == length(ω_list) == length(T_list)) || throw(DimensionMismatch("The number of fields, κs, ωs and Ts must be the same."))

    dims = N_trunc == size(H,1) ? H.dims : [N_trunc]
    E, U = eigen(H)
    U = QuantumObject(U, dims=H.dims)

    H_d = QuantumObject(dense_to_sparse((U' * H * U).data[1:N_trunc, 1:N_trunc], tol), dims=dims)

    Ω = dense_to_sparse((E' .- E)[1:N_trunc,1:N_trunc], tol)
    Ω = triu(Ω, 1)

    L = liouvillian(H_d)

    for i in eachindex(fields)
        X_op = dense_to_sparse((U' * fields[i] * U).data[1:N_trunc, 1:N_trunc], tol)

        # # Sp₀ = QuantumObject( (Ω ./ ω_list[i]) .* droptol!(sparse(triu(X_op, 1)), tol) )
        # Sp₀ = QuantumObject( droptol!(sparse(triu(X_op, 1)), tol) )
        # Sp₁ = QuantumObject( n_th.(Ω, T_list[i]) .* Sp₀.data )
        # Sp2₂ = QuantumObject( (1 .+ n_th.(Ω, T_list[i])) .* Sp₀.data )

        # Nikki Master Equation
        N_th = n_th.(Ω, T_list[i])
        Sp₀ = QuantumObject( triu(X_op, 1), dims=dims )
        Sp₁ = QuantumObject( (Ω ./ ω_list[i]) .* N_th .* Sp₀.data, dims=dims )
        Sp2₂ = QuantumObject( (Ω ./ ω_list[i]) .* (1 .+ N_th) .* Sp₀.data, dims=dims )
        S0 = QuantumObject( spdiagm(diag(X_op)), dims=dims )

        L += κ_list[i] / 2 * ( sprepost(Sp₁', Sp₀) + sprepost(Sp₀', Sp₁) - spre(Sp₀ * Sp₁') - spost(Sp₁ * Sp₀') )
        L += κ_list[i] / 2 * ( sprepost(Sp2₂, Sp₀') + sprepost(Sp₀, Sp2₂') - spre(Sp₀' * Sp2₂) - spost(Sp2₂' * Sp₀) )
        L += κ_list[i] * T_list[i] / (4 * ω_list[i]) * ( 4 * sprepost(S0, S0) - 2 * spre(S0 * S0) - 2 * spost(S0 * S0) )
    end

    return E, U, L
end

function _liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    Lₚ::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    Lₘ::QuantumObject{<:AbstractArray{T3},SuperOperatorQuantumObject},
    ω::Real, solver::LiouvillianDirectSolver; n_max::Int=4) where {T1,T2,T3}

    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data

    L_0_d = sparse_to_dense(L_0)
    L_p_d = sparse_to_dense(L_p)
    L_m_d = sparse_to_dense(L_m)

    n_i = n_max
    S, T = -(L_0_d - 1im * n_i * ω * I) \ L_p_d, -(L_0_d + 1im * n_i * ω * I) \ L_m_d

    for n_i in n_max-1:-1:1
        S, T = -(L_0_d - 1im * n_i * ω * I + L_m_d * S) \ L_p_d, -(L_0_d + 1im * n_i * ω * I + L_p_d * T) \ L_m_d
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

function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T1},OpType1},
    H_p::QuantumObject{<:AbstractArray{T2},OpType2},
    H_m::QuantumObject{<:AbstractArray{T3},OpType3},
    ω::Real; n_max::Int=4, lf_solver::Type{LSolver}=LiouvillianDirectSolver,
    ss_solver::Type{SSSolver}=SteadyStateDirectSolver) where {T1,T2,T3,OpType1<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType2<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    OpType3<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    LSolver<:LiouvillianSolver,SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end