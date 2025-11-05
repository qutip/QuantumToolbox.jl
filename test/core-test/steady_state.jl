@testitem "Steady State" begin
    N = 10
    a = destroy(N)
    a_d = a'
    H = a_d * a + 0.1 * (a + a_d)
    Ht = QobjEvo(H, (p, t) -> 1) # for test throw
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
    rho_me = sol_me.states[end]

    solver = SteadyStateODESolver()
    ρ_ss = steadystate(H, c_ops, solver = solver)
    @test tracedist(rho_me, ρ_ss) < 1e-4

    solver = SteadyStateDirectSolver()
    ρ_ss = steadystate(H, c_ops, solver = solver)
    @test tracedist(rho_me, ρ_ss) < 1e-4
    @test_throws ArgumentError steadystate(Ht, c_ops, solver = solver)

    solver = SteadyStateLinearSolver()
    ρ_ss = steadystate(H, c_ops, solver = solver)
    @test tracedist(rho_me, ρ_ss) < 1e-4
    @test_throws ArgumentError steadystate(Ht, c_ops, solver = solver)

    solver = SteadyStateEigenSolver()
    ρ_ss = steadystate(H, c_ops, solver = solver)
    @test tracedist(rho_me, ρ_ss) < 1e-4
    @test_throws ArgumentError steadystate(Ht, c_ops, solver = solver)

    @testset "Type Inference (steadystate)" begin
        L = liouvillian(H, c_ops)

        solver = SteadyStateODESolver(tmax = t_l[end])
        @inferred steadystate(H, c_ops, solver = solver)
        @inferred steadystate(L, solver = solver)

        solver = SteadyStateDirectSolver()
        @inferred steadystate(H, c_ops, solver = solver)
        @inferred steadystate(L, solver = solver)

        solver = SteadyStateLinearSolver()
        @inferred steadystate(H, c_ops, solver = solver)
        @inferred steadystate(L, solver = solver)

        solver = SteadyStateEigenSolver()
        @inferred steadystate(H, c_ops, solver = solver)
        @inferred steadystate(L, solver = solver)

        @inferred steadystate(H, c_ops)
        @inferred steadystate(L)
    end

    H = a_d * a
    H_t = 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 100 * 2π, 1000)

    coeff(p, t) = sin(t)
    H_td = (H, (H_t, coeff))

    sol_me = mesolve(H_td, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
    ρ_ss1 = steadystate_fourier(H, -1im * 0.5 * H_t, 1im * 0.5 * H_t, 1, c_ops, solver = SteadyStateLinearSolver())[1]
    ρ_ss2 =
        steadystate_fourier(H, -1im * 0.5 * H_t, 1im * 0.5 * H_t, 1, c_ops, solver = SSFloquetEffectiveLiouvillian())

    @test abs(sum(sol_me.expect[1, (end-100):end]) / 101 - expect(e_ops[1], ρ_ss1)) < 1e-3
    @test abs(sum(sol_me.expect[1, (end-100):end]) / 101 - expect(e_ops[1], ρ_ss2)) < 1e-3

    @testset "Type Inference (steadystate_fourier)" begin
        @inferred steadystate_fourier(
            H,
            -1im * 0.5 * H_t,
            1im * 0.5 * H_t,
            1,
            c_ops,
            solver = SteadyStateLinearSolver(),
        )
        @inferred steadystate_fourier(
            H,
            -1im * 0.5 * H_t,
            1im * 0.5 * H_t,
            1,
            c_ops,
            solver = SSFloquetEffectiveLiouvillian(),
        )
    end
end
