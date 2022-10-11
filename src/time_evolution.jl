function LindbladJumpCallback(savebefore=false,saveafter=false)
    function LindbladJumpCondition(u, t, integrator)
        norm(u)^2 - integrator.p[2]
    end

    function LindbladJumpAffect!(integrator)
        ψ = integrator.u
        c_ops = integrator.p[1]
        
        if length(c_ops) == 1
            integrator.u = normalize(c_ops[1] * ψ)
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
                prob += real(expect(c_op' * c_op, ψ)) / dp
                if prob >= r2
                    collaps_idx = i
                    break
                end
            end
            integrator.u = normalize(c_ops[collaps_idx] * ψ)
        end
        integrator.p = [c_ops, rand()]
    end

    ContinuousCallback(LindbladJumpCondition, LindbladJumpAffect!, save_positions = (savebefore,saveafter))
end

"""
    function mcsolve(H::AbstractArray, ψ0, t_l, c_ops; 
        e_ops = [], 
        n_traj::Int = 1,
        batch_size::Int = n_traj % 10 == 0 ? round(Int, n_traj / 10) : n_traj,
        alg = AutoVern7(KenCarp4(autodiff=false)),
        ensemble_method = EnsembleThreads(), 
        update_function = (t)->0*I,
        progress = true,
        kwargs...)

Time evolution of an open quantum system using quantum trajectories.
"""
function mcsolve(H::AbstractArray, ψ0, t_l, c_ops;
    e_ops = [], 
    n_traj::Int = 1,
    batch_size::Int = n_traj % 10 == 0 ? round(Int, n_traj / 10) : n_traj,
    alg = AutoVern7(KenCarp4(autodiff=false)),
    ensemble_method = EnsembleThreads(), 
    update_function = (t)->0*I,
    progress = true,
    kwargs...)

    tspan = (t_l[1], t_l[end])

    H_eff = H
    for c_op in c_ops
        H_eff += - 0.5im * c_op' * c_op
    end

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
        res = hcat(map(i->map(op->expect(op, normalize(sol.u[i])), e_ops), eachindex(t_l))...)
        (res, false)
    end
    function reduction(u,batch,I)
        tmp = sum(cat(batch..., dims = 3), dims = 3)
        length(u) == 0 && return tmp, false
        cat(u, tmp, dims = 3), false
    end

    dudt!(du,u,p,t) = mul!(du, -1im * (H_eff + update_function(t)), u)

    cb1 = LindbladJumpCallback()
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2)

    p = [c_ops, rand()]
    prob = ODEProblem(dudt!, ψ0, tspan, p, callback = cb; kwargs...)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, reduction=reduction)
    sol = solve(ensemble_prob, alg, ensemble_method, trajectories=n_traj, batch_size=batch_size, saveat = t_l)

    put!(channel, false)

    length(e_ops) == 0 && return sol

    e_ops_expect = sum(sol.u, dims = 3) ./ n_traj

    return sol, e_ops_expect
end

"""
    mesolve(H::AbstractArray, ψ0, t_l, c_ops; 
        e_ops = [], 
        alg = Vern7(), 
        update_function = nothing, 
        progress = true,
        kwargs...)

Time evolution of an open quantum system using master equation.
"""
function mesolve(H::AbstractArray, ψ0, t_l, c_ops; 
    e_ops = [], 
    alg = Vern7(), 
    update_function = nothing, 
    progress = true,
    kwargs...)

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ρ0 = reshape(ψ0 * ψ0', length(ψ0)^2)

    L = liouvillian(H, c_ops)
    function L_t(t)
        update_function === nothing && return 0 * I
        size(update_function(0)) == size(H) && (ft = update_function(t); return -1im * (spre(ft) - spost(ft)))
        return update_function(t)
    end

    saved_values = SavedValues(Float64, Vector{Float64})
    function save_func(u, t, integrator)
        next!(progr)
        map(op->expect(op, reshape(u, size(H)...)), e_ops)
    end
    cb1 = SavingCallback(save_func, saved_values, saveat = t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2)

    if typeof(alg) <: LinearExponential
        !(update_function === nothing) && error("The Liouvillian must to be time independent when using LinearExponential algorith.")
        A = DiffEqArrayOperator(L)
        prob = ODEProblem(A, 
        ρ0, tspan; kwargs...)
        sol = solve(prob, alg, dt = (t_l[2] - t_l[1]), callback = cb)
    else
        dudt!(du,u,p,t) = mul!(du, L + L_t(t), u)
        prob = ODEProblem(dudt!, 
        ρ0, tspan; kwargs...)
        sol = solve(prob, alg, callback = cb)
    end

    length(e_ops) == 0 && return sol

    return sol, hcat(saved_values.saveval...)
end

"""
    sesolve(H::AbstractArray, ψ0, t_l; 
        e_ops = [], 
        alg = Vern7(), 
        update_function = nothing, 
        progress = true,
        kwargs...)

Time evolution of a closed quantum system using Schrödinger equation.
"""
function sesolve(H::AbstractArray, ψ0, t_l; 
    e_ops = [], 
    alg = Vern7(), 
    update_function = nothing, 
    progress = true,
    kwargs...)

    H_t(t) = update_function === nothing ? 0*I : update_function(t)

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    saved_values = SavedValues(Float64, Vector{Float64}) 
    function save_func(u, t, integrator)
        next!(progr)
        map(op->expect(op, u), e_ops)
    end
    cb1 = SavingCallback(save_func, saved_values, saveat = t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2)

    if typeof(alg) <: LinearExponential
        !(update_function === nothing) && error("The Hamiltonian must to be time independent when using LinearExponential algorithm.")
        A = DiffEqArrayOperator(-1im * H)
        prob = ODEProblem(A, ψ0, tspan, callback = cb; kwargs...)
        sol = solve(prob, alg, dt = (t_l[2] - t_l[1]))
    else
        dudt!(du,u,p,t) = mul!(du, -1im * (H + H_t(t)), u)
        prob = ODEProblem(dudt!, ψ0, tspan, callback = cb; kwargs...)
        sol = solve(prob, alg)
    end

    length(e_ops) == 0 && return sol

    return sol, hcat(saved_values.saveval...)
end

function liouvillian(H::AbstractArray, c_ops)
    L = -1im * (spre(H) - spost(H))
    for c_op in c_ops
        L += lindblad_dissipator(c_op)
    end
    L
end

liouvillian(H::AbstractArray) = -1im * (spre(H) - spost(H))

function liouvillian_floquet(L_0::AbstractArray, L_p::AbstractArray, L_m::AbstractArray, w_l::Real; n_max::Int = 4)
    N_size = size(L_0, 1)
    Id = eye(N_size)
    S = T = spzeros(eltype(L_0), N_size, N_size)

    L_p_d = Matrix(L_p)
    L_m_d = Matrix(L_m)

    for n_i in n_max:-1:1
        S, T = - ( L_0 - 1im * n_i * w_l * Id + L_m * S ) \ L_p_d, - ( L_0 + 1im * n_i * w_l * Id + L_p * T ) \ L_m_d
    end

    return droptol!(sparse(L_0 + L_m * S + L_p * T), 1e-8)
end

function steadystate(L::AbstractArray)
    L_tmp = copy(L)
    N_size = floor(Int, √(size(L_tmp, 1)))
    weight = sum( abs.(L_tmp) ) / length(L_tmp)
    v0 = zeros(ComplexF64, N_size^2)
    v0[1] = weight
    
    L_tmp[1, [N_size * (i - 1) + i for i in 1:N_size]] .+= weight
    
    rho_ss_vec = L_tmp \ v0
    rho_ss = sparse(reshape(rho_ss_vec, N_size, N_size))
    return rho_ss
end

function steadystate(H::AbstractArray, c_ops::Vector)
    L = liouvillian(H, c_ops)

    return steadystate(L)
end

function steadystate_floquet(H_0::AbstractArray, c_ops::Vector, H_p::AbstractArray, H_m::AbstractArray, w_l::Real; n_max::Int = 4)
    L_0 = liouvillian(H_0, c_ops)
    L_p = liouvillian(H_p)
    L_m = liouvillian(H_m)

    steadystate(liouvillian_floquet(L_0, L_p, L_m, w_l, n_max = n_max))
end