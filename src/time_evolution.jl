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
            for i in eachindex(c_ops)
                c_op = c_ops[i]
                dp += real(psi' * (c_op' * c_op) * psi)
            end
            prob = 0
            for i in eachindex(c_ops)
                c_op = c_ops[i]
                prob += real(psi' * (c_op' * c_op) * psi) / dp
                if prob >= r2[]
                    collaps_idx = i
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
function mcsolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], n_traj = 1, ensemble_method = EnsembleSerial(), update_function = nothing, abstol = 1e-7, reltol = 1e-5)
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

    if update_function === nothing
        A = DiffEqArrayOperator(-1im * H_eff)
    else
        function update_func(A, u, p, t)
            copyto!(A, -1im * (H_eff + update_function(t)))
        end
        
        A = DiffEqArrayOperator(-1im * H_eff, update_func = update_func)
    end
    prob = ODEProblem(A, psi0, tspan, p)
    cb1 = LindbladJumpCallback()
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2)

    prob = ODEProblem{true}(A, psi0, tspan, p, callback = cb, abstol = abstol, reltol = reltol)
    ensemble_prob = EnsembleProblem(prob)
    sol = solve(ensemble_prob, alg_hints=[:nonstiff], ensemble_method, trajectories=n_traj, saveat = t_l)

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))
    end

    Threads.@threads for i in 1:length(t_l)
        for idx in eachindex(sol)
            normalize!(sol[idx].u[i])
        end
        for j in eachindex(e_ops)
            e_ops_expect[i, j] = mean([expect(e_ops[j], sol[idx2].u[i]) for idx2 in eachindex(sol)])
        end
    end

    if length(e_ops) == 0
        return sol
    else
        return sol, e_ops_expect
    end
end

function mesolve(H::AbstractArray, psi0, t_l, c_ops; e_ops = [], krylovdim = 10, abstol = 1e-7, reltol = 1e-5)
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
    prob = ODEProblem(A, rho0_vec, tspan, abstol = abstol, reltol = reltol)
    sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (N_t - 1))

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

        Threads.@threads for i in 1:length(t_l)
            for j in eachindex(e_ops)
                e_ops_expect[i, j] = expect(e_ops[j], reshape(sol.u[i], size(H, 1), size(H, 1)))
            end
        end

        return sol, e_ops_expect
    else
        return sol
    end
end

function sesolve(H::AbstractArray, psi0, t_l; e_ops = [], krylovdim = 10, abstol = 1e-7, reltol = 1e-5)
    N_t = length(t_l)
    ti, tf = t_l[1], t_l[end]
    tspan = (ti, tf)

    A = DiffEqArrayOperator(-1im * H)
    prob = ODEProblem(A, psi0, tspan, abstol = abstol, reltol = reltol)
    sol = solve(prob, LinearExponential(krylov=:adaptive, m = krylovdim), dt = (tf - ti) / (N_t - 1))

    if length(e_ops) != 0
        e_ops_expect = zeros(Float64, length(t_l), length(e_ops))

        Threads.@threads for i in 1:length(t_l)
            for j in eachindex(e_ops)
                e_ops_expect[i, j] = expect(e_ops[j], sol.u[i])
            end
        end

        return sol, e_ops_expect
    else
        return sol
    end
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