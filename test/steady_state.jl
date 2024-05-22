@testset "Steady State" begin
    N = 10
    a = destroy(N)
    a_d = a'
    H = a_d * a + 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = false)
    ρ_ss = steadystate(H, c_ops)
    @test abs(sol_me.expect[1, end] - expect(e_ops[1], ρ_ss)) < 1e-3

    H = a_d * a
    H_t = 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    H_t_f = TimeDependentOperatorSum([(t, p) -> sin(t)], [liouvillian(H_t)])
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, H_t = H_t_f, alg = Vern7(), progress_bar = false)
    ρ_ss = steadystate_floquet(H, c_ops, -1im * 0.5 * H_t, 1im * 0.5 * H_t, 1)
    @test abs(sum(sol_me.expect[1, end-100:end]) / 101 - expect(e_ops[1], ρ_ss)) < 1e-2
end
