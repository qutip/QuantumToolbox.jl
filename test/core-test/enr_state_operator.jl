@testitem "Excitation number restricted state space" begin
    @testset "mesolve, steadystate, and eigenstates" begin
        ε = 2π
        ωc = 2π
        g = 0.1ωc
        γ = 0.01ωc
        tlist = range(0, 20, 100)
        N_cut = 2

        # normal mesolve and steadystate
        sz = sigmaz() ⊗ qeye(N_cut)
        sm = sigmam() ⊗ qeye(N_cut)
        a = qeye(2) ⊗ destroy(N_cut)
        H_JC = 0.5ε * sz + ωc * a' * a + g * (sm' * a + a' * sm)
        c_ops_JC = (√γ * a,)
        ψ0_JC = basis(2, 0) ⊗ fock(N_cut, 0)
        sol_JC = mesolve(H_JC, ψ0_JC, tlist, c_ops_JC; e_ops = [sz], progress_bar = Val(false))
        ρ_ss_JC = steadystate(H_JC, c_ops_JC)

        # ENR mesolve and steadystate
        N_exc = 1
        dims = (2, N_cut)
        sm_enr, a_enr = enr_destroy(dims, N_exc)
        sz_enr = 2 * sm_enr' * sm_enr - 1
        ψ0_enr = enr_fock(dims, N_exc, [1, 0])
        H_enr = ε * sm_enr' * sm_enr + ωc * a_enr' * a_enr + g * (sm_enr' * a_enr + a_enr' * sm_enr)
        c_ops_enr = (√γ * a_enr,)
        sol_enr = mesolve(H_enr, ψ0_enr, tlist, c_ops_enr; e_ops = [sz_enr], progress_bar = Val(false))
        ρ_ss_enr = steadystate(H_enr, c_ops_enr)
        
        # check mesolve result
        @test all(isapprox.(sol_JC.expect, sol_enr.expect, atol = 1e-4))

        # check steadystate result
        @test expect(sz, ρ_ss_JC) ≈ expect(sz_enr, ρ_ss_enr) atol=1e-4

        # check eigenstates
        λ, v = eigenstates(H_enr)
        @test all([H_enr * v[k] ≈ λ[k] * v[k] for k in eachindex(λ)])

        qeye(3) ⊗ H_enr
    end
end
