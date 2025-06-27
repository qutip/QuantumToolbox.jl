@testitem "Excitation number restricted state space" begin
    @testset "mesolve" begin
        ε = 2π
        ωc = 2π
        g = 0.1ωc
        γ = 0.01ωc
        tlist = range(0, 20, 100)
        N_cut = 2

        sz = sigmaz() ⊗ qeye(N_cut)
        sm = sigmam() ⊗ qeye(N_cut)
        a = qeye(2) ⊗ destroy(N_cut)

        H_JC = 0.5ε * sz + ωc * a' * a + g * (sm' * a + a' * sm)
        ψ0 = basis(2, 0) ⊗ fock(N_cut, 0)
        sol_JC = mesolve(H_JC, ψ0, tlist, [√γ * a]; e_ops = [sz], progress_bar = Val(false))

        N_exc = 1
        dims = [2, N_cut]
        sm_enr, a_enr = enr_destroy(dims, N_exc)
        sz_enr = 2 * sm_enr' * sm_enr - 1
        ψ0_enr = enr_fock(dims, N_exc, [1, 0])
        H_enr = ε * sm_enr' * sm_enr + ωc * a_enr' * a_enr + g * (sm_enr' * a_enr + a_enr' * sm_enr)
        sol_enr = mesolve(H_enr, ψ0_enr, tlist, [√γ * a_enr]; e_ops = [sz_enr], progress_bar = Val(false))

        @test all(isapprox.(sol_JC.expect, sol_enr.expect, atol = 1e-4))
    end
end
