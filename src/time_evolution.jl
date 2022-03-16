function LindbladJumpCallback()
    r1 = Ref{Float64}(rand(Float64))
    r2 = Ref{Float64}(rand(Float64))
    condition = (u, t, integrator) -> real(u' * u) - r1[]

    affect! = function (integrator)
        psi = copy(integrator.u)
        c_ops = copy(integrator.p)
        collaps_idx = 1
        if length(c_ops) == 1
            integrator.u = normalize(c_ops[1] * psi)
        else
            r2[] = rand()
            dp = 0
            for i in 1:length(c_ops)
                c_op = c_ops[i]
                dp += real(psi' * (c_op' * c_op) * psi)
            end
            prob = 0
            for i in 1:length(c_ops)
                c_op = c_ops[i]
                prob += real(psi' * (c_op' * c_op) * psi) / dp
                # println(prob)
                if prob >= r2[]
                    collaps_idx = i
                    # println(collaps_idx)
                    break
                end
            end
            integrator.u = normalize(c_ops[collaps_idx] * psi)
        end
        
        r1[] = rand()
    end

    initialize = function (c, u, t, integrator)
        r1[] = rand()
        u_modified!(integrator, false)
    end

    return ContinuousCallback(condition, affect!, initialize = initialize, save_positions = (false,false))
end

"""
    mcsolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], n_traj = 100, ensemble_method = EnsembleSerial())

Time evolution of an open quantum system using quantum trajectories.
"""
function mcsolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], n_traj = 100, ensemble_method = EnsembleSerial())
    N_c_ops = length(c_ops)
    N_t = length(t_l)
    ti, tf = t_l[1], t_l[end]
    tspan = (ti, tf)
    p = c_ops

    H_eff = H
    for i in 1:N_c_ops
        c_op = c_ops[i]
        H_eff += - 0.5im * c_op' * c_op
    end
    
    A = DiffEqArrayOperator(-1im * H_eff)
    prob = ODEProblem(A, psi0, tspan, p)
    cb1 = LindbladJumpCallback()
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2)

    prob = ODEProblem{true}(A, psi0, tspan, p, callback = cb, abstol = 1e-7, reltol = 1e-5)
    ensemble_prob = EnsembleProblem(prob)
    # Tsit5() alg_hints=[:nonstiff] lsoda()
    sol = solve(ensemble_prob, alg_hints=[:nonstiff], ensemble_method, trajectories=n_traj, saveat = t_l)
    # sol = solve(ensemble_prob, LinearExponential(krylov=:adaptive), 
    #     EnsembleThreads(), trajectories=n_traj, dt = (tf - ti) / N_t)

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))
    end

    Threads.@threads for i in 1:length(t_l)
        for idx in 1:length(sol)
            normalize!(sol[idx].u[i])
        end
        for j in 1:length(e_ops)
            e_ops_expect[i, j] = mean([expect(e_ops[j], sol[idx2].u[i]) for idx2 in 1:length(sol)])
        end
    end

    # Threads.@threads for idx in 1:length(sol)
    #     for i in 1:length(sol[idx].t)
    #         normalize!(sol[idx].u[i])
    #     end
    # end

    if length(e_ops) == 0
        return sol
    else
        return sol, e_ops_expect
    end
end

function mesolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], krylovdim = 10)
    N_c_ops = length(c_ops)
    N_t = length(t_l)
    ti, tf = t_l[1], t_l[end]
    tspan = (ti, tf)
    p = c_ops

    L = -1im * (spre(H) - spost(H))
    for i in 1:N_c_ops
        L += lindblad_dissipator( c_ops[i] )
    end

    rho0_vec = reshape(psi0 * psi0', length(psi0)^2)

    A = DiffEqArrayOperator(L)
    prob = ODEProblem(A, rho0_vec, tspan, abstol = 1e-7, reltol = 1e-5)
    sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (N_t - 1))

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

        Threads.@threads for i in 1:length(t_l)
            for j in 1:length(e_ops)
                if ishermitian(e_ops[j])
                    e_ops_expect[i, j] = real(tr(e_ops[j] * reshape(sol.u[i], size(H, 1), size(H, 1))))
                else
                    e_ops_expect[i, j] = tr(e_ops[j] * reshape(sol.u[i], size(H, 1), size(H, 1)))
                end
            end
        end

        return sol, e_ops_expect
    else
        return sol
    end
end

function sesolve(H::AbstractArray, psi0, t_l; e_ops = [], krylovdim = 10)
    N_t = length(t_l)
    ti, tf = t_l[1], t_l[end]
    tspan = (ti, tf)

    A = DiffEqArrayOperator(-1im * H)
    prob = ODEProblem(A, psi0, tspan, abstol = 1e-7, reltol = 1e-5)
    sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (N_t - 1))

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

        Threads.@threads for i in 1:length(t_l)
            for j in 1:length(e_ops)
                e_ops_expect[i, j] = expect(e_ops[j], sol.u[i])
            end
        end

        return sol, e_ops_expect
    else
        return sol
    end
end