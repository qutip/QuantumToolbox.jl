@testset "Dynamical Shifted Fock" begin
    F = 3
    Δ = 0.25
    κ = 1
    U = 0.01

    tlist = LinRange(0, 25, 300)

    # Single cavity case
    N0 = 100
    a0 = destroy(N0)
    H0 = Δ * a0' * a0 + F * (a0 + a0') + U * a0'^2 * a0^2
    c_ops0 = [√(κ) * a0]

    α0 = 1.5
    ρ0 = coherent(N0, α0)
    sol0 = mesolve(H0, ρ0, tlist, c_ops0, e_ops = [a0' * a0, a0], progress_bar = false)

    N = 5
    a = destroy(N)
    function H_dsf(op_list, p)
        Δ = p.Δ
        F = p.F
        U = p.U
        a = op_list[1]
        return Δ * a' * a + F * (a + a') + U * a'^2 * a^2
    end
    function c_ops_dsf(op_list, p)
        κ = p.κ
        a = op_list[1]
        return [√κ * a]
    end
    function e_ops_dsf(op_list, p)
        a = op_list[1]
        return [a' * a, a]
    end
    op_list = [a]
    ψ0 = fock(N, 0)
    α0_l = [α0]
    dsf_params = (Δ = Δ, F = F, κ = κ, U = U)

    sol_dsf_me =
        dsf_mesolve(H_dsf, ψ0, tlist, c_ops_dsf, op_list, α0_l, dsf_params, e_ops = e_ops_dsf, progress_bar = false)
    sol_dsf_mc = dsf_mcsolve(
        H_dsf,
        ψ0,
        tlist,
        c_ops_dsf,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf,
        progress_bar = false,
        n_traj = 500,
    )
    val_ss = abs2(sol0.expect[1, end])
    @test sum(abs2.(sol0.expect[1, :] .- sol_dsf_me.expect[1, :])) / (val_ss * length(tlist)) < 0.1
    @test sum(abs2.(sol0.expect[1, :] .- sol_dsf_mc.expect[1, :])) / (val_ss * length(tlist)) < 0.1

    # Two cavities case
    F = 2
    Δ = 0.25
    κ = 1
    U = 0.01
    J = 0.5
    tlist = LinRange(0, 15, 300)

    N0 = 20
    a10 = kron(destroy(N0), qeye(N0))
    a20 = kron(qeye(N0), destroy(N0))
    H0 =
        Δ * a10' * a10 +
        Δ * a20' * a20 +
        U * a10'^2 * a10^2 +
        U * a20'^2 * a20^2 +
        F * (a10 + a10') +
        J * (a10' * a20 + a10 * a20')
    c_ops0 = [√κ * a10, √κ * a20]

    ρ0 = kron(coherent(N0, α0), coherent(N0, α0))
    sol0 = mesolve(H0, ρ0, tlist, c_ops0, e_ops = [a10' * a10, a20' * a20], progress_bar = false)

    N = 5
    a1 = kron(destroy(N), qeye(N))
    a2 = kron(qeye(N), destroy(N))
    function H_dsf2(op_list, p)
        Δ = p.Δ
        F = p.F
        U = p.U
        J = p.J
        a1, a2 = op_list
        return Δ * a1' * a1 +
               Δ * a2' * a2 +
               U * a1'^2 * a1^2 +
               U * a2'^2 * a2^2 +
               F * (a1 + a1') +
               J * (a1' * a2 + a1 * a2')
    end
    function c_ops_dsf2(op_list, p)
        κ = p.κ
        a1, a2 = op_list
        return [√κ * a1, √κ * a2]
    end
    function e_ops_dsf2(op_list, p)
        a1, a2 = op_list
        return [a1' * a1, a2' * a2]
    end
    op_list = [a1, a2]
    ψ0 = kron(fock(N, 0), fock(N, 0))
    α0_l = [α0, α0]
    dsf_params = (Δ = Δ, F = F, κ = κ, U = U, J = J)

    sol_dsf_me =
        dsf_mesolve(H_dsf2, ψ0, tlist, c_ops_dsf2, op_list, α0_l, dsf_params, e_ops = e_ops_dsf2, progress_bar = false)
    sol_dsf_mc = dsf_mcsolve(
        H_dsf2,
        ψ0,
        tlist,
        c_ops_dsf2,
        op_list,
        α0_l,
        dsf_params,
        e_ops = e_ops_dsf2,
        progress_bar = false,
        n_traj = 500,
    )

    val_ss = abs2(sol0.expect[1, end])
    @test sum(abs2.(sol0.expect[1, :] .- sol_dsf_me.expect[1, :])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[1, :] .- sol_dsf_mc.expect[1, :])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[2, :] .- sol_dsf_me.expect[2, :])) / (val_ss * length(tlist)) < 0.6
    @test sum(abs2.(sol0.expect[2, :] .- sol_dsf_mc.expect[2, :])) / (val_ss * length(tlist)) < 0.6
end
