using QuPhys
using Test


@testset "QuPhys.jl" begin
    # Write your tests here.
end

@testset "Time Evolution and partial trace" begin
    N = 10

    a_d = kron(create(N), eye(2))
    a = a_d'
    sm = kron(eye(N), sigmam())
    sp = sm'
    sx = sm + sp
    sy = 1im * (sm - sp)
    sz = sp * sm - sm * sp
    η = 0.01
    H = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx
    psi0 = kron(fock(N, 0), fock(2, 0))
    t_l = LinRange(0, 1000, 1000)
    e_ops = [a_d * a]
    sol, expect_se = sesolve(H, psi0, t_l, e_ops = e_ops, progress = false)
    @test sum(abs.(expect_se[1, :] .- sin.(η * t_l).^2)) / length(t_l) < 0.1
    sol, expect_se = sesolve(H, psi0, t_l, e_ops = e_ops, alg = Vern7(), progress = false)
    @test sum(abs.(expect_se[1, :] .- sin.(η * t_l).^2)) / length(t_l) < 0.1

    a = destroy(N)
    a_d = a'
    H = a_d * a
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 100, 1000)
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress = false);
    sol_mc, expect_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = e_ops, progress = false);
    @test sum(abs.(expect_mc .- expect_me)) / length(t_l) < 0.1
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, alg = Vern7(), progress = false);
    @test sum(abs.(expect_mc .- expect_me)) / length(t_l) < 0.1

    sp1 = kron(sigmap(), eye(2))
    sm1 = sp1'
    sx1 = sm1 + sp1
    sy1 = 1im * (sm1 - sp1)
    sz1 = sp1 * sm1 - sm1 * sp1
    sp2 = kron(eye(2), sigmap())
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
    t_l = LinRange(0, 10 / (γ1 + γ2), 1000)
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = [sp1 * sm1, sp2 * sm2], progress = false);
    sol_mc, expect_mc = mcsolve(H, psi0, t_l, c_ops, n_traj = 500, e_ops = [sp1 * sm1, sp2 * sm2], progress = false);
    @test sum(abs.(expect_mc[1:2, :] .- expect_me[1:2, :])) / length(t_l) < 0.1

    ## partial trace
    @test expect(sp1 * sm1, reshape(sol_me(4 / (γ1 + γ2)), 4, 4)) ≈ expect(sigmap() * sigmam(), ptrace(reshape(sol_me(4 / (γ1 + γ2)), 4, 4), [1], (2, 2)))
end

@testset "Steadystate" begin
    N = 10

    a = destroy(N)
    a_d = a'
    H = a_d * a + 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress = false);
    ρ_ss = steadystate(H, c_ops)
    @test abs(expect_me[1, end] - expect(e_ops[1], ρ_ss)) < 1e-3

    H = a_d * a
    H_t = 0.1 * (a + a_d)
    c_ops = [sqrt(0.1) * a]
    e_ops = [a_d * a]
    psi0 = fock(N, 3)
    t_l = LinRange(0, 200, 1000)
    sol_me, expect_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, update_function = (t)->sin(t)*H_t, progress = false);
    ρ_ss = steadystate_floquet(H, c_ops, -1im * 0.5 * H_t, 1im * 0.5 * H_t, 1)
    @test abs(sum(expect_me[1, end-100:end])/101 - expect(e_ops[1], ρ_ss)) < 1e-2
end

@testset "Entanglement" begin
    g = fock(2, 1)
    e = fock(2, 0)
    state = normalize(kron(g, e) + kron(e, g))
    rho = state * state'
    @test entanglement(state, [1], (2, 2)) / log(2) ≈ 1
end

@testset "Wigner" begin
    α = 0.5 + 0.8im
    ψ = coherent(30, α)
    xvec = LinRange(-3, 3, 300)
    yvec = LinRange(-3, 3, 300)

    wig = wigner(ψ, xvec, yvec)

    X, Y = meshgrid(xvec, yvec)
    wig_tmp1 = gaussian(xvec / √2, real(α), 1 / 2)
    wig_tmp2 = gaussian(yvec / √2, imag(α), 1 / 2)
    wig2 = maximum(wig) * reshape(kron(wig_tmp1, wig_tmp2), 300, 300)

    @test sqrt(sum(abs.(wig2 .- wig)) / length(wig)) < 0.1
end