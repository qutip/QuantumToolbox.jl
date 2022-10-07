function LindbladJumpCallback()
    function LindbladJumpCondition(u, t, integrator)
        norm(u)^2 < integrator.p[2]
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
            for i in eachindex(c_ops)
                c_op = c_ops[i]
                dp += real(ψ' * (c_op' * c_op) * ψ)
            end
            prob = 0
            for i in eachindex(c_ops)
                c_op = c_ops[i]
                prob += real(ψ' * (c_op' * c_op) * ψ) / dp
                if prob >= r2
                    collaps_idx = i
                    break
                end
            end
            integrator.u = normalize(c_ops[collaps_idx] * ψ)
        end
        integrator.p = [c_ops, rand()]
    end

    DiscreteCallback(LindbladJumpCondition, LindbladJumpAffect!, save_positions = (false,false))
end

"""
    mcsolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], n_traj = 1, ensemble_method = EnsembleSerial(), update_function = (t)->0*I, abstol = 1e-7, reltol = 1e-5)

Time evolution of an open quantum system using quantum trajectories.
"""
function mcsolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], n_traj = 1, ensemble_method = EnsembleSerial(), update_function = (t)->0*I, kwargs...)
    tspan = (t_l[1], t_l[end])

    H_eff = H
    for c_op in c_ops
        H_eff += - 0.5im * c_op' * c_op
    end

    dudt!(du,u,p,t) = mul!(du, -1im * (H_eff + update_function(t)), u)

    cb1 = LindbladJumpCallback()

    function prob_func(prob,i,repeat)
        remake(prob,p=[prob.p[1], rand()])
    end

    p = [c_ops, rand()]
    prob = ODEProblem(dudt!, psi0, tspan, p, kwargs...)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    sol = solve(ensemble_prob, Vern7(), callback = cb1, ensemble_method, trajectories=n_traj, saveat = t_l)

    length(e_ops) == 0 && return sol

    e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

    for i in 1:length(t_l)
        for j in eachindex(e_ops)
            e_ops_expect[i, j] = sum([expect(e_ops[j], normalize!(sol[idx].u[i])) for idx in eachindex(sol)]) / length(sol)
        end
    end

    return sol, e_ops_expect
end

"""
    mesolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], method = "Vern7", update_function = nothing, krylovdim = 30, kwargs...)

Time evolution of an open quantum system using master equation.
"""
function mesolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], method = "Vern7", update_function = nothing, krylovdim = 30, kwargs...)
    tspan = (t_l[1], t_l[end])

    L = -1im * (spre(H) - spost(H))
    for c_op in c_ops
        L += lindblad_dissipator( c_op )
    end
    function L_t(t)
        update_function === nothing && return 0 * I
        size(update_function(0)) == size(H) && (ft = update_function(t); return -1im * (spre(ft) - spost(ft)))
        return update_function(t)
    end

    rho0_vec = reshape(psi0 * psi0', length(psi0)^2)

    if method == "LinearExponential"
        !(update_function === nothing) && error("The Liouvillian must to be time independent when using LinearExponential method.")
        A = DiffEqArrayOperator(L)
        prob = ODEProblem(A, rho0_vec, tspan, kwargs...)
        sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (length(t_l) - 1))
    elseif method == "Vern7"
        dudt!(du,u,p,t) = mul!(du, L + L_t(t), u)
        prob = ODEProblem(dudt!, rho0_vec, tspan, kwargs...)
        sol = solve(prob, Vern7(), saveat = t_l)
    end

    length(e_ops) == 0 && return sol

    e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

    for i in eachindex(t_l)
        for j in eachindex(e_ops)
            e_ops_expect[i, j] = expect(e_ops[j], reshape(sol.u[i], size(H)...))
        end
    end

    return sol, e_ops_expect
end

"""
    sesolve(H::AbstractArray, psi0, t_l; e_ops = [], method = "Vern7", update_function = nothing, krylovdim = 10, kwargs...)

Time evolution of a closed quantum system using Schrödinger equation.
"""
function sesolve(H::AbstractArray, psi0, t_l; e_ops = [], method = "Vern7", update_function = nothing, krylovdim = 10, kwargs...)
    tspan = (t_l[1], t_l[end])

    function H_t(t)
        update_function === nothing && return 0 * I
        return update_function(t)
    end

    if method == "LinearExponential"
        !(update_function === nothing) && error("The Hamiltonian must to be time independent when using LinearExponential method.")
        A = DiffEqArrayOperator(-1im * H)
        prob = ODEProblem(A, psi0, tspan, kwargs...)
        sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (length(t_l) - 1))
    elseif method == "Vern7"
        dudt!(du,u,p,t) = mul!(du, -1im * (H + H_t(t)), u)
        prob = ODEProblem(dudt!, psi0, tspan, kwargs...)
        sol = solve(prob, Vern7(), saveat = t_l)
    end

    length(e_ops) == 0 && return sol

    e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

    for i in 1:length(t_l)
        for j in eachindex(e_ops)
            e_ops_expect[i, j] = expect(e_ops[j], sol.u[i])
        end
    end

    return sol, e_ops_expect
end

function liouvillian(H::AbstractArray, c_ops)
    L = -1im * (spre(H) - spost(H))
    for c_op in c_ops
        L += lindblad_dissipator(c_op)
    end
    L
end

function liouvillian_floquet(L_0::AbstractArray, L_p::AbstractArray, L_m::AbstractArray, w_l; n_max = 4)
    N_size = size(L_0, 1)
    Id = eye(N_size)
    S = T = spzeros(ComplexF64, N_size, N_size)

    L_p_d = Matrix(L_p)
    L_m_d = Matrix(L_m)

    for n_i in n_max:-1:1
        S = - ( L_0 - 1im * n_i * w_l * Id + L_m * S ) \ L_p_d
        T = - ( L_0 + 1im * n_i * w_l * Id + L_p * T ) \ L_m_d
    end

    return L_0 + L_m * S + L_p * T
end

function steadystate(L::AbstractArray)
    L_tmp = copy(L)
    N_size = floor(Int, √(size(L_tmp, 1)))
    weight = mean( abs.(L_tmp) )
    v0 = zeros(ComplexF64, N_size^2)
    v0[1] = weight
    
    L_tmp[1, [N_size * (i - 1) + i for i in 1:N_size]] .+= weight
    
    rho_ss_vec = L_tmp \ v0
    rho_ss = sparse(reshape(rho_ss_vec, N_size, N_size))
    return rho_ss
end

function steadystate(H::AbstractArray, c_ops::Vector)
    L = -1im * ( spre(H) - spost(H) )
    for op in c_ops
        L += lindblad_dissipator(op)
    end

    return steadystate(L)
end