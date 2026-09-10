@testitem "Excitation number restricted state space (evolution)" begin
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
    @test all(isapprox.(sol_JC.expect, sol_enr.expect, atol = 1.0e-4))

    # check steadystate result
    @test expect(sz, ρ_ss_JC) ≈ expect(sz_enr, ρ_ss_enr) atol = 1.0e-4

    # check eigenstates
    λ, v = eigenstates(H_enr)
    @test all([H_enr * v[k] ≈ λ[k] * v[k] for k in eachindex(λ)])
end

@testitem "Excitation number restricted state space (weighted evolution)" begin
    # Parametric down-conversion conserves the weighted number n_a + n_b + 2 n_c, but not
    # the total number, so it requires `excitation_weights`. Compare the full Fock-space
    # dynamics with the (much smaller) ENR-space dynamics.
    Δa, Δb, Δc = 1.0, 1.3, 2.1
    g = 0.5
    γ = 0.1
    tlist = range(0, 10, 50)

    # full Fock space (signal, idler, pump)
    Na = Nb = Nc = 2
    a = destroy(Na) ⊗ qeye(Nb) ⊗ qeye(Nc)
    b = qeye(Na) ⊗ destroy(Nb) ⊗ qeye(Nc)
    c = qeye(Na) ⊗ qeye(Nb) ⊗ destroy(Nc)
    Mf = a' * b' * c
    H = Δa * a' * a + Δb * b' * b + Δc * c' * c + g * (Mf + Mf')
    ψ0 = fock(Na, 0) ⊗ fock(Nb, 0) ⊗ fock(Nc, 1)  # one pump excitation, N = 2
    c_ops = (√γ * a, √γ * b)
    sol_full = mesolve(H, ψ0, tlist, c_ops; e_ops = [a' * a, c' * c], progress_bar = Val(false))

    # ENR space with conserved weighted number n_a + n_b + 2 * n_c <= 2
    dims = (2, 2, 2)
    n_exc = 2
    weights = (1, 1, 2)
    s_enr = EnrSpace(dims, n_exc; excitation_weights = weights)
    ae, be, ce = enr_destroy(s_enr)
    Me = ae' * be' * ce  # annihilation on the right => intermediate states stay in the ENR space
    H_enr = Δa * ae' * ae + Δb * be' * be + Δc * ce' * ce + g * (Me + Me')
    ψ0_enr = enr_fock(s_enr, [0, 0, 1])
    c_ops_enr = (√γ * ae, √γ * be)
    sol_enr = mesolve(H_enr, ψ0_enr, tlist, c_ops_enr; e_ops = [ae' * ae, ce' * ce], progress_bar = Val(false))

    @test size(H_enr.data, 1) == 5  # ENR truncation (vs 8 for the full space)
    @test all(isapprox.(sol_full.expect, sol_enr.expect, atol = 1.0e-5))
end
