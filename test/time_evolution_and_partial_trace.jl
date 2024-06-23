@testset "Time Evolution and Partial Trace" begin
    # sesolve
    N = 10
    a_d = kron(create(N), qeye(2))
    a = a_d'
    sx = kron(qeye(N), sigmax())
    sy = tensor(qeye(N), sigmay())
    sz = qeye(N) ⊗ sigmaz()
    η = 0.01
    H = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx
    psi0 = kron(fock(N, 0), fock(2, 0))
    t_l = LinRange(0, 1000, 1000)
    e_ops = [a_d * a]
    sol = sesolve(H, psi0, t_l, e_ops = e_ops, alg = Vern7(), progress_bar = false)
    sol2 = sesolve(H, psi0, t_l, progress_bar = false)
    sol3 = sesolve(H, psi0, t_l, e_ops = e_ops, saveat = t_l, progress_bar = false)
    sol_string = sprint((t, s) -> show(t, "text/plain", s), sol)
    @test sum(abs.(sol.expect[1, :] .- sin.(η * t_l) .^ 2)) / length(t_l) < 0.1
    @test ptrace(sol.states[end], 1) ≈ ptrace(ket2dm(sol.states[end]), 1)
    @test ptrace(sol.states[end]', 1) ≈ ptrace(sol.states[end], 1)
    @test length(sol.states) == 1
    @test size(sol.expect) == (length(e_ops), length(t_l))
    @test length(sol2.states) == length(t_l)
    @test size(sol2.expect) == (0, length(t_l))
    @test length(sol3.states) == length(t_l)
    @test size(sol3.expect) == (length(e_ops), length(t_l))
    @test sol_string ==
          "Solution of time evolution\n" *
          "(return code: $(sol.retcode))\n" *
          "--------------------------\n" *
          "num_states = $(length(sol.states))\n" *
          "num_expect = $(size(sol.expect, 1))\n" *
          "ODE alg.: $(sol.alg)\n" *
          "abstol = $(sol.abstol)\n" *
          "reltol = $(sol.reltol)\n"

    # mesolve and mcsolve
    a = destroy(N)
    a_d = a'
    H = a_d * a
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = basis(N, 3)
    t_l = LinRange(0, 100, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, alg = Vern7(), progress_bar = false)
    sol_me2 = mesolve(H, psi0, t_l, c_ops, progress_bar = false)
    sol_me3 = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, saveat = t_l, progress_bar = false)
    sol_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = e_ops, progress_bar = false)
    sol_me_string = sprint((t, s) -> show(t, "text/plain", s), sol_me)
    sol_mc_string = sprint((t, s) -> show(t, "text/plain", s), sol_mc)
    @test sum(abs.(sol_mc.expect .- sol_me.expect)) / length(t_l) < 0.1
    @test length(sol_me.states) == 1
    @test size(sol_me.expect) == (length(e_ops), length(t_l))
    @test length(sol_me2.states) == length(t_l)
    @test size(sol_me2.expect) == (0, length(t_l))
    @test length(sol_me3.states) == length(t_l)
    @test size(sol_me3.expect) == (length(e_ops), length(t_l))
    @test sol_me_string ==
          "Solution of time evolution\n" *
          "(return code: $(sol_me.retcode))\n" *
          "--------------------------\n" *
          "num_states = $(length(sol_me.states))\n" *
          "num_expect = $(size(sol_me.expect, 1))\n" *
          "ODE alg.: $(sol_me.alg)\n" *
          "abstol = $(sol_me.abstol)\n" *
          "reltol = $(sol_me.reltol)\n"
    @test sol_mc_string ==
          "Solution of quantum trajectories\n" *
          "(converged: $(sol_mc.converged))\n" *
          "--------------------------------\n" *
          "num_trajectories = $(sol_mc.n_traj)\n" *
          "num_states = $(length(sol_mc.states[1]))\n" *
          "num_expect = $(size(sol_mc.expect, 1))\n" *
          "ODE alg.: $(sol_mc.alg)\n" *
          "abstol = $(sol_mc.abstol)\n" *
          "reltol = $(sol_mc.reltol)\n"

    # exceptions
    psi_wrong = basis(N - 1, 3)
    @test_throws DimensionMismatch sesolve(H, psi_wrong, t_l)
    @test_throws DimensionMismatch mesolve(H, psi_wrong, t_l)
    @test_throws DimensionMismatch mcsolve(H, psi_wrong, t_l)
    @test_throws ArgumentError sesolve(H, psi0, t_l, save_idxs = [1, 2])
    @test_throws ArgumentError mesolve(H, psi0, t_l, save_idxs = [1, 2])
    @test_throws ArgumentError mcsolve(H, psi0, t_l, save_idxs = [1, 2])

    sp1 = kron(sigmap(), qeye(2))
    sm1 = sp1'
    sx1 = sm1 + sp1
    sy1 = 1im * (sm1 - sp1)
    sz1 = sp1 * sm1 - sm1 * sp1
    sp2 = kron(qeye(2), sigmap())
    sm2 = sp2'
    sx2 = sm2 + sp2
    sy2 = 1im * (sm2 - sp2)
    sz2 = sp2 * sm2 - sm2 * sp2
    ωq1, ωq2 = 1, 1
    γ1, γ2 = 0.05, 0.1
    H = 0.5 * ωq1 * sz1 + 0.5 * ωq2 * sz2
    c_ops = [sqrt(γ1) * sm1, sqrt(γ2) * sm2]
    psi0_1 = normalize(fock(2, 0) + fock(2, 1))
    psi0_2 = normalize(fock(2, 0) + fock(2, 1))
    psi0 = kron(psi0_1, psi0_2)
    t_l = LinRange(0, 20 / γ1, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = [sp1 * sm1, sp2 * sm2], progress_bar = false)
    sol_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = [sp1 * sm1, sp2 * sm2], progress_bar = false)
    @test sum(abs.(sol_mc.expect[1:2, :] .- sol_me.expect[1:2, :])) / length(t_l) < 0.1
    @test expect(sp1 * sm1, sol_me.states[end]) ≈ expect(sigmap() * sigmam(), ptrace(sol_me.states[end], 1))
end
