abstract type LiouvillianSolver end
abstract type LiouvillianDirectSolver <: LiouvillianSolver end

abstract type SteadyStateSolver end
abstract type SteadyStateDirectSolver <: SteadyStateSolver end

struct TimeEvolutionSol
    times::AbstractVector
    states::AbstractArray
    expect::AbstractArray
end

function LindbladJumpCallback(savebefore=false,saveafter=false)
    function LindbladJumpCondition(u, t, integrator)
        norm(u)^2 - integrator.p[2]
    end

    function LindbladJumpAffect!(integrator)
        ψ = integrator.u
        c_ops = integrator.p[1]
        
        if length(c_ops) == 1
            integrator.u = normalize!(c_ops[1] * ψ)
        else
            collaps_idx = 1
            r2 = rand()
            dp = 0
            @inbounds for i in eachindex(c_ops)
                c_op = c_ops[i]
                dp += real(ψ' * (c_op' * c_op) * ψ)
            end
            prob = 0
            @inbounds for i in eachindex(c_ops)
                c_op = c_ops[i]
                prob += real(ψ' * c_op' * c_op * ψ) / dp
                if prob >= r2
                    collaps_idx = i
                    break
                end
            end
            integrator.u = normalize!(c_ops[collaps_idx] * ψ)
        end
        integrator.p = [c_ops, rand()]
    end

    ContinuousCallback(LindbladJumpCondition, nothing, LindbladJumpAffect!, interp_points=0, save_positions = (savebefore,saveafter))
end

"""
    function mcsolve(H::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, 
            ψ0::QuantumObject{<:AbstractArray{T}, KetQuantumObject}, 
            t_l::AbstractVector, c_ops::AbstractVector;
            e_ops::AbstractVector = [], 
            n_traj::Int = 1,
            batch_size::Int = min(10, n_traj),
            alg = AutoVern7(KenCarp4(autodiff=false)),
            ensemble_method = EnsembleThreads(), 
            H_t = nothing,
            progress = true,
            callbacks = [],
            kwargs...)

Time evolution of an open quantum system using quantum trajectories.
"""
function mcsolve(H::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, 
            ψ0::QuantumObject{<:AbstractArray{T}, KetQuantumObject}, 
            t_l::AbstractVector, c_ops::AbstractVector;
            e_ops::AbstractVector = [], 
            n_traj::Int = 1,
            batch_size::Int = min(10, n_traj),
            alg = AutoVern7(KenCarp4(autodiff=false)),
            ensemble_method = EnsembleThreads(), 
            H_t = nothing,
            progress = true,
            callbacks = [],
            kwargs...) where {T}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims

    tspan = (t_l[1], t_l[end])
    e_ops_len = length(e_ops)

    H_eff = H
    for c_op in c_ops
        H_eff += - 0.5im * c_op' * c_op
    end
    H_eff = H_eff.data
    ψ0 = ψ0.data
    c_ops = [op.data for op in c_ops]

    progr = Progress(n_traj, showspeed=true, enabled=progress)
    channel = RemoteChannel(()->Channel{Bool}(), 1)
    @async while take!(channel)
        next!(progr)
    end

    function prob_func(prob,i,repeat)
        remake(prob,p=[prob.p[1], rand()])
    end
    function output_func(sol,i)
        put!(channel, true)
        if e_ops_len == 0
            # res = hcat(sol.u...)
            res = [QuantumObject(ϕ, dims=Hdims) for ϕ in sol.u]
        else
            res = hcat(map(i->map(op->expect(op, QuantumObject(normalize!(sol.u[i]), dims=Hdims)), e_ops), eachindex(t_l))...)
        end
        (res, false)
    end
    function reduction(u,batch,I)
        if e_ops_len == 0
            tmp = hcat(batch...)
            length(u) == 0 && return tmp, false
            res = hcat(u, tmp)
        else
            tmp = sum(cat(batch..., dims = 3), dims = 3)
            length(u) == 0 && return tmp, false
            res = sum(cat(u, tmp, dims = 3), dims = 3)
        end
        return res, false
    end

    is_time_dependent = !(H_t === nothing)
    if is_time_dependent
        dudt! = (du,u,p,t) -> mul!(du, -1im * (H_eff + H_t(t).data), u)
    else
        dudt! = (du,u,p,t) -> mul!(du, -1im * H_eff, u)
    end

    cb1 = LindbladJumpCallback()
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, callbacks...)

    p = [c_ops, rand()]
    prob = ODEProblem(dudt!, ψ0, tspan, p, callback = cb; kwargs...)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, reduction=reduction)
    sol = solve(ensemble_prob, alg, ensemble_method, trajectories=n_traj, 
            batch_size=batch_size, saveat = t_l)

    put!(channel, false)

    e_ops_len == 0 && return TimeEvolutionSol(t_l, sol.u, [])

    e_ops_expect = sum(sol.u, dims = 3) ./ n_traj

    return TimeEvolutionSol(t_l, [], e_ops_expect)
end

"""
    function mesolve(H::QuantumObject{<:AbstractArray{T}, HOpType}, 
            ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
            t_l::AbstractVector, c_ops::AbstractVector; 
            e_ops::AbstractVector = [], 
            alg = LinearExponential(krylov=:simple), 
            H_t = nothing, 
            params::AbstractVector = [],
            progress = true,
            callbacks = [],
            kwargs...)

Time evolution of an open quantum system using master equation.
"""
function mesolve(H::QuantumObject{<:AbstractArray{T}, HOpType}, 
            ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
            t_l::AbstractVector, c_ops::AbstractVector; 
            e_ops::AbstractVector = [], 
            alg = LinearExponential(krylov=:adaptive, m=15),
            H_t = nothing,
            params::AbstractVector = [],
            progress = true,
            callbacks = [],
            kwargs...) where {T,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
                                StateOpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims
    Hsize = prod(Hdims)

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ψ0_data = ψ0.data
    if isket(ψ0)
        ρ0 = reshape(ψ0_data * ψ0_data', length(ψ0)^2)
    elseif isbra(ψ0)
        ρ0 = reshape(ψ0_data' * ψ0_data, length(ψ0)^2)
    else
        ρ0 = reshape(ψ0_data, length(ψ0))
    end

    L = liouvillian(H, c_ops).data

    is_time_dependent = !(H_t === nothing)

    saved_values = SavedValues(Float64, Vector{ComplexF64})
    function save_func(u, t, integrator)
        next!(progr)
        map(op->expect(op, QuantumObject(reshape(u, Hsize, Hsize), OperatorQuantumObject, Hdims)), e_ops)
    end
    cb1 = SavingCallback(save_func, saved_values, saveat = t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, callbacks...)

    if typeof(alg) <: LinearExponential
        is_time_dependent && error("The Liouvillian must to be time independent when using LinearExponential algorith.")
        A = DiffEqArrayOperator(L)
        prob = ODEProblem(A, ρ0, tspan, params; kwargs...)
        sol = solve(prob, alg, dt = (t_l[2] - t_l[1]), callback = cb)
    else
        if !is_time_dependent
            dudt! = (du,u,p,t) -> mul!(du, L, u)
        else
            if H_t(0.0).type <: OperatorQuantumObject
                @warn string("To speed up the calculation, it is always better to define ",
                    "the time-dependent part as a SuperOperator, and not as an Operator.") maxlog=1
                dudt! = (du,u,p,t) -> mul!(du, L + liouvillian(H_t(t)).data, u)
            else
                dudt! = (du,u,p,t) -> mul!(du, L + H_t(t).data, u)
            end
        end
        prob = ODEProblem(dudt!, ρ0, tspan, params; kwargs...)
        sol = solve(prob, alg, callback = cb)
    end

    ρt_len = isqrt(length(sol.u[1]))
    if ρt_len > 1
        ρt = [QuantumObject(sparse(reshape(ϕ, ρt_len, ρt_len)), dims=Hdims) for ϕ in sol.u]
    else
        ρt = []
    end
    length(e_ops) == 0 && return TimeEvolutionSol(sol.t, ρt, [])

    return TimeEvolutionSol(sol.t, ρt, hcat(saved_values.saveval...))
end

"""
    function sesolve(H::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject},
                ψ0::QuantumObject{<:AbstractArray{T}, KetQuantumObject},
                t_l::AbstractVector;  
                e_ops::AbstractVector = [], 
                alg = LinearExponential(), 
                H_t = nothing, 
                params::AbstractVector = [],
                progress = true,
                callbacks = [],
                kwargs...)

Time evolution of a closed quantum system using Schrödinger equation.
"""
function sesolve(H::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject},
            ψ0::QuantumObject{<:AbstractArray{T}, KetQuantumObject},
            t_l::AbstractVector;  
            e_ops::AbstractVector = [], 
            alg = LinearExponential(), 
            H_t = nothing, 
            params::AbstractVector = [],
            progress = true,
            callbacks = [],
            kwargs...) where {T}

    H.dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H.dims

    tspan = (t_l[1], t_l[end])

    H = H.data
    ψ0 = ψ0.data

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    is_time_dependent = !(H_t === nothing)

    saved_values = SavedValues(Float64, Vector{Float64}) 
    function save_func(u, t, integrator)
        next!(progr)
        map(op->expect(op, QuantumObject(normalize!(u), dims=Hdims)), e_ops)
    end
    cb1 = SavingCallback(save_func, saved_values, saveat = t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, callbacks...)

    if typeof(alg) <: LinearExponential
        is_time_dependent && error("The Hamiltonian must to be time independent when using LinearExponential algorithm.")
        A = DiffEqArrayOperator(-1im * H)
        prob = ODEProblem(A, ψ0, tspan, params; kwargs...)
        sol = solve(prob, alg, dt = (t_l[2] - t_l[1]), callback = cb)
    else
        if !is_time_dependent
            dudt! = (du,u,p,t) -> mul!(du, -1im * H, u)
        else
            dudt! = (du,u,p,t) -> mul!(du, -1im * (H + H_t(t).data), u)
        end
        prob = ODEProblem(dudt!, ψ0, tspan, params; kwargs...)
        sol = solve(prob, alg, callback = cb)
    end

    ψt_len = isqrt(length(sol.u[1]))
    if ψt_len == prod(Hdims)
        ψt = [QuantumObject(ϕ, dims=Hdims) for ϕ in sol.u]
    else
        ψt = []
    end
    length(e_ops) == 0 && return TimeEvolutionSol(sol.t, ψt, [])

    return TimeEvolutionSol(sol.t, ψt, hcat(saved_values.saveval...))
end

### DYNAMICAL FOCK DIMENSION ###
function _reduce_dims(QO::QuantumObject{<:AbstractArray{T}, OpType}, sel::AbstractVector, reduce::AbstractVector) where {T,OpType<:OperatorQuantumObject}
    rd = QO.dims
    reduce_l = zero(rd)
    reduce_l[sel] .= reduce
    rd_new = rd .- reduce_l

    ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
    ρmat2 = zeros(eltype(ρmat), repeat(rd_new, 2)...)
    copyto!(ρmat2, view(ρmat, repeat([1:n for n in rd_new], 2)...))
    ρmat = reshape(PermutedDimsArray(ρmat2, length(size(ρmat2)):-1:1), prod(rd_new), prod(rd_new))

    QuantumObject(ρmat, OperatorQuantumObject, rd_new)
end

function _increase_dims(QO::QuantumObject{<:AbstractArray{T}, OpType}, sel::AbstractVector, increase::AbstractVector) where {T,OpType<:OperatorQuantumObject}
    rd = QO.dims
    incr_l = zero(rd)
    incr_l[sel] .= increase
    rd_new = rd .+ incr_l

    ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
    ρmat2 = zeros(eltype(ρmat), repeat(rd_new, 2)...)
    copyto!(view(ρmat2, repeat([1:n for n in rd], 2)...), ρmat)
    ρmat = reshape(PermutedDimsArray(ρmat2, length(size(ρmat2)):-1:1), prod(rd_new), prod(rd_new))

    QuantumObject(ρmat, OperatorQuantumObject, rd_new)
end

function _DFDIncreaseReduceCondition(u, t, integrator)
    internal_params = integrator.p[1]
    dim_list = internal_params["dim_list"]
    maxdims = internal_params["maxdims"]
    tol_list = internal_params["tol_list"]
    increase_list = internal_params["increase_list"]
    reduce_list = internal_params["reduce_list"]
    pillow_list = internal_params["pillow_list"]
    # Here I use copy(u) otherwise I have an error.
    ρt = QuantumObject(reshape(copy(u), prod(dim_list), prod(dim_list)), OperatorQuantumObject, dim_list)

    condition = false
    @inbounds for i in 1:length(dim_list)
        maxdim_i = maxdims[i]
        dim_i = dim_list[i]
        if dim_i < maxdim_i && dim_i > 2 && maxdim_i != 0
            ρi = ptrace(ρt, [i]).data
            pillow_i = min( max(round(Int, 0.02 * dim_i), 1), 20 )
            pillow_list[i] = pillow_i
            @views res = sum(abs.(ρi[diagind(ρi)[end-pillow_i:end]])) / pillow_i
            if res > tol_list[i]
                push!(increase_list, i)
                condition = true
            elseif res < tol_list[i]*1e-2 && dim_i > 3
                push!(reduce_list, i)
                condition = true
            end
        end
    end
    condition
end

function _DFDIncreaseReduceAffect!(integrator)
    internal_params = integrator.p[1]
    H = internal_params["H"]
    c_ops = internal_params["c_ops"]
    dim_list = internal_params["dim_list"]
    increase_list = internal_params["increase_list"]
    reduce_list = internal_params["reduce_list"]
    pillow_list = internal_params["pillow_list"]
    # Here I use copy(integrator.u) otherwise I have an error.
    ρt  = QuantumObject(reshape(copy(integrator.u), prod(dim_list), prod(dim_list)), OperatorQuantumObject, dim_list)

    @views pillow_increase = pillow_list[increase_list]
    @views pillow_reduce = pillow_list[reduce_list]
    if length(increase_list) > 0
        ρt = _increase_dims(ρt, increase_list, pillow_increase)
        dim_list[increase_list] .+= pillow_increase
    end
    if length(reduce_list) > 0
        ρt = _reduce_dims(ρt, reduce_list, pillow_reduce)
        dim_list[reduce_list] .-= pillow_reduce
    end

    deleteat!(increase_list, 1:length(increase_list))
    deleteat!(reduce_list, 1:length(reduce_list))

    params = integrator.p
    L = liouvillian(H(dim_list), c_ops(dim_list)).data
    newsize = size(L, 1)
    resize!(integrator, newsize)
    params[1]["L"] = L
    integrator.u = reshape(ρt.data, newsize)
end

"""
    function dfd_mesolve(H::Function, ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
        t_l::AbstractVector, c_ops::Function, e_ops::Function, maxdims::AbstractVector;
        tol_list::AbstractVector = [],
        alg = Vern7(),
        progress::Bool = true,
        callbacks::AbstractVector = [],
        kwargs...)

Time evolution of an open quantum system using master equation, dynamically changing the dimension of the Hilbert subspaces.
"""
function dfd_mesolve(H::Function, ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
    t_l::AbstractVector, c_ops::Function, e_ops::Function, maxdims::AbstractVector;
    tol_list::AbstractVector = [],
    alg = Vern7(),
    progress::Bool = true,
    callbacks::AbstractVector = [],
    kwargs...) where {T,StateOpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    dim_list = copy(ψ0.dims)
    length(dim_list) != length(maxdims) && throw(DimensionMismatch("'dim_list' and 'maxdims' do not have the same dimension."))
    H(dim_list).dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ψ0_data = ψ0.data
    if isket(ψ0)
        ρ0 = reshape(ψ0_data * ψ0_data', length(ψ0)^2)
    elseif isbra(ψ0)
        ρ0 = reshape(ψ0_data' * ψ0_data, length(ψ0)^2)
    else
        ρ0 = reshape(ψ0_data, length(ψ0))
    end
    
    L = liouvillian(H(dim_list), c_ops(dim_list)).data
    
    length(tol_list) == 0 && (tol_list = [1e-8 for d in dim_list])
    reduce_list = Int16[]
    increase_list = Int16[]
    pillow_list = ones(Int16, length(dim_list))

    saved_values = SavedValues(Float64, Vector{ComplexF64})
    function save_func(u, t, integrator)
        internal_params = integrator.p[1]
        e_ops = internal_params["e_ops"]
        dim_list = internal_params["dim_list"]
        next!(progr)
        # Here I use copy(u) otherwise I have an error.
        ρt = QuantumObject(reshape(copy(u), prod(dim_list), prod(dim_list)), OperatorQuantumObject, dim_list)
        map(op->expect(op, ρt), e_ops(dim_list))
    end
    cb1 = SavingCallback(save_func, saved_values, saveat = t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb_increasereduce = DiscreteCallback(_DFDIncreaseReduceCondition, _DFDIncreaseReduceAffect!, save_positions=(false,false))
    cb = CallbackSet(cb1, cb2, cb_increasereduce, callbacks...)
    
    params = [Dict("L" => L, "H" => H, "c_ops" => c_ops, "e_ops" => e_ops, 
                   "dim_list" => dim_list, "maxdims" => maxdims, "tol_list" => tol_list, 
                   "reduce_list" => reduce_list, "increase_list" => increase_list, "pillow_list" => pillow_list)]
    
    dudt! = (du,u,p,t) -> mul!(du, p[1]["L"], u)
    prob = ODEProblem(dudt!, ρ0, tspan, params; kwargs...)
    sol = solve(prob, alg, callback = cb)

    ρt = [QuantumObject(reshape(ϕ, isqrt(length(ϕ)), :)) for ϕ in sol.u]

    TimeEvolutionSol(sol.t, ρt, hcat(saved_values.saveval...))
end


### LIOUVILLIAN AND STEADYSTATE ###
function liouvillian(H::QuantumObject{<:AbstractArray{T}, OpType}, 
            c_ops::AbstractVector) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}

    L = isoper(H) ? -1im * (spre(H) - spost(H)) : H
    for c_op in c_ops
        if isoper(c_op)
            L += lindblad_dissipator(c_op)
        elseif issuper(c_op)
            L += c_op
        end
    end
    L
end


liouvillian(H::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = 
liouvillian(H, [])

function liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            Lₚ::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            Lₘ::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            ω::Real; n_max::Int = 4, solver::Type{LSolver} = LiouvillianDirectSolver) where {T1,LSolver<:LiouvillianSolver}

    ((L₀.dims == Lₚ.dims) && (L₀.dims == Lₘ.dims)) || throw(ErrorException("The operators are not of the same Hilbert dimension."))

    _liouvillian_floquet(L₀, Lₚ, Lₘ, ω, solver, n_max=n_max)
end

function _liouvillian_floquet(L₀::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            Lₚ::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            Lₘ::QuantumObject{<:AbstractArray{T1}, SuperOperatorQuantumObject}, 
            ω::Real, solver::Type{LiouvillianDirectSolver}; n_max::Int = 4) where {T1}

    L_0 = L₀.data
    L_p = Lₚ.data
    L_m = Lₘ.data
    S = T = spzeros(T1, size(L_0)...)

    L_p_d = Matrix(L_p)
    L_m_d = Matrix(L_m)

    for n_i in n_max:-1:1
        S, T = - ( L_0 - 1im * n_i * ω * I + L_m_d * S ) \ L_p_d, - ( L_0 + 1im * n_i * ω * I + L_p_d * T ) \ L_m_d
    end

    QuantumObject(droptol!(sparse(L_0 + L_m * S + L_p * T), 1e-12), SuperOperatorQuantumObject, L₀.dims)
end

function steadystate(L::QuantumObject{<:AbstractArray{T}, SuperOperatorQuantumObject}; 
            solver::Type{SSSolver} = SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    _steadystate(L, solver)
end

function steadystate(H::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, c_ops::Vector, 
            solver::Type{SSSolver} = SteadyStateDirectSolver) where {T,SSSolver<:SteadyStateSolver}

    L = liouvillian(H, c_ops)
    steadystate(L, solver=solver)
end

function _steadystate(L::QuantumObject{<:AbstractArray{T}, SuperOperatorQuantumObject}, 
            solver::Type{SteadyStateDirectSolver}) where {T}

    L_tmp = copy(L.data)
    N = prod(L.dims) # floor(Int, √(size(L_tmp, 1)))
    weight = sum( abs.(L_tmp) ) / length(L_tmp)
    v0 = zeros(ComplexF64, N^2)
    v0[1] = weight

    L_tmp[1, [N * (i - 1) + i for i in 1:N]] .+= weight

    rho_ss_vec = L_tmp \ v0
    rho_ss = droptol!(sparse(reshape(rho_ss_vec, N, N)), 1e-12)
    QuantumObject(rho_ss, OperatorQuantumObject, L.dims)
end

function steadystate_floquet(H_0::QuantumObject{<:AbstractArray{T}, OpType1}, 
            c_ops::Vector, H_p::QuantumObject{<:AbstractArray{T}, OpType2}, 
            H_m::QuantumObject{<:AbstractArray{T}, OpType3}, 
            ω::Real; n_max::Int = 4, lf_solver::Type{LSolver} = LiouvillianDirectSolver,
            ss_solver::Type{SSSolver} = SteadyStateDirectSolver) where {T,OpType1<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                        OpType2<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                        OpType3<:Union{OperatorQuantumObject, SuperOperatorQuantumObject},
                                                                        LSolver<:LiouvillianSolver, SSSolver<:SteadyStateSolver}

    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, ω, n_max=n_max, solver=lf_solver), solver=ss_solver)
end