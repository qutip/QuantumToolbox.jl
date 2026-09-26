using Test
using QuantumToolbox

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "Time-dependent Hamiltonian" begin
    # Get parameters from TESetup to simplify the code
    ωd = TESetup.ωd
    F = TESetup.F
    a = TESetup.a
    σz = TESetup.σz
    H = TESetup.H
    H_td = TESetup.H_td
    H_td2 = TESetup.H_td2
    L_td = TESetup.L_td
    ψ0 = TESetup.ψ0
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    p = TESetup.p
    rng = TESetup.rng

    # ssesolve is slow to be run on CI. It is not removed from the test because it may be useful for testing in more powerful machines.

    # Time Evolution in the drive frame

    H_dr_fr = H - ωd * a' * a - ωd * σz / 2 + F * (a + a')

    tlist = range(0, 10 / TESetup.γ, 1000)

    sol_se = sesolve(H_dr_fr, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
    sol_me = mesolve(H_dr_fr, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
    sol_mc = mcsolve(H_dr_fr, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)
    # sol_sse = ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)

    # Time Evolution in the lab frame

    sol_se_td = sesolve(H_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_me_td = mesolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_mc_td = mcsolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)
    # sol_sse_td = ssesolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)

    @test sol_se.expect ≈ sol_se_td.expect atol = 1.0e-6 * length(tlist)
    @test sol_me.expect ≈ sol_me_td.expect atol = 1.0e-6 * length(tlist)
    @test sol_mc.expect ≈ sol_mc_td.expect atol = 1.0e-2 * length(tlist)
    # @test sol_sse.expect ≈ sol_sse_td.expect atol = 1e-2 * length(tlist)

    sol_se_td2 = sesolve(H_td2, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_me_td2 = mesolve(L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_mc_td2 = mcsolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)
    # sol_sse_td2 =
    # ssesolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)

    @test sol_se.expect ≈ sol_se_td2.expect atol = 1.0e-6 * length(tlist)
    @test sol_me.expect ≈ sol_me_td2.expect atol = 1.0e-6 * length(tlist)
    @test sol_mc.expect ≈ sol_mc_td2.expect atol = 1.0e-2 * length(tlist)
    # @test sol_sse.expect ≈ sol_sse_td2.expect atol = 1e-2 * length(tlist)

    # test assume_hermitian = Val(false)
    L_td_non_herm = liouvillian(H_td2, c_ops, assume_hermitian = Val(false))
    sol_me_td3 = mesolve(L_td_non_herm, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
    @test sol_me.expect ≈ sol_me_td3.expect atol = 1.0e-6 * length(tlist)
end

@testset "Time-dependent collapse operators" begin
    rng = TESetup.rng

    # ---- Analytic benchmark: a single qubit with a time-dependent decay rate ----
    # For C(t) = √γ(t) σ⁻, the excited-state population obeys dP/dt = -γ(t) P, hence
    # P(t) = exp(-∫₀ᵗ γ(t')dt'). With γ(t) = γ₀(1 + ½cos(ωt)) (strictly positive),
    # ∫₀ᵗ γ = γ₀(t + ½ sin(ωt)/ω).
    σ = sigmam()
    pq = (γ0 = 0.4, ω = 0.5)
    fq(p, t) = sqrt(p.γ0 * (1 + 0.5 * cos(p.ω * t)))
    c_op_q = QobjEvo(σ, fq)
    Hq = 0 * qeye(2)
    ψe = basis(2, 0)  # excited state
    Pe = σ' * σ       # excited-state projector
    tlist_q = range(0, 10, 100)
    analytic = @. exp(-pq.γ0 * (tlist_q + 0.5 * sin(pq.ω * tlist_q) / pq.ω))

    sol_me_q = mesolve(Hq, ψe, tlist_q, (c_op_q,), e_ops = (Pe,), params = pq, progress_bar = Val(false))
    @test real.(sol_me_q.expect[1, :]) ≈ analytic atol = 1.0e-4 # mesolve is essentially exact here

    sol_mc_q =
        mcsolve(Hq, ψe, tlist_q, (c_op_q,), e_ops = (Pe,), params = pq, ntraj = 1000, rng = rng, progress_bar = Val(false))
    @test real.(sol_mc_q.expect[1, :]) ≈ analytic atol = 1.0e-2 * length(tlist_q)

    # ---- mesolve ↔ mcsolve consistency on a richer system (decay + heating) ----
    N = TESetup.N
    a = TESetup.a
    σm = TESetup.σm
    H = TESetup.H
    e_ops = TESetup.e_ops
    γ = TESetup.γ
    nth = TESetup.nth

    ψ0 = kron(fock(N, 1), fock(2, 1)) # start with excitation so the dynamics are nontrivial
    tlist = range(0, 10 / γ, 100)

    # Time-dependent decay (`f_dn`) and heating (`f_up`) rates; `p` is passed as `params`.
    # The heating channels keep the total jump rate strictly positive (no dark state),
    # mirroring the constant `c_ops` of the test setup.
    p = (γ = γ, nth = nth, ω = 0.5)
    f_dn(p, t) = sqrt(p.γ * (1 + p.nth) * (1 + 0.5 * sin(p.ω * t)^2))
    f_up(p, t) = sqrt(p.γ * p.nth * (1 + 0.5 * sin(p.ω * t)^2))
    c_ops_td = (QobjEvo(a, f_dn), QobjEvo(a', f_up), QobjEvo(σm, f_dn), QobjEvo(σm', f_up))

    sol_me = mesolve(H, ψ0, tlist, c_ops_td, e_ops = e_ops, params = p, progress_bar = Val(false))
    sol_mc = mcsolve(H, ψ0, tlist, c_ops_td, e_ops = e_ops, params = p, ntraj = 500, rng = rng, progress_bar = Val(false))
    @test sol_me.expect ≈ sol_mc.expect atol = 1.0e-2 * length(tlist)

    # A mix of constant and time-dependent collapse operators must also work
    c_ops_mixed = (QobjEvo(a, f_dn), sqrt(γ * nth) * a', sqrt(γ * (1 + nth)) * σm, sqrt(γ * nth) * σm')
    sol_me_mix = mesolve(H, ψ0, tlist, c_ops_mixed, e_ops = e_ops, params = p, progress_bar = Val(false))
    sol_mc_mix =
        mcsolve(H, ψ0, tlist, c_ops_mixed, e_ops = e_ops, params = p, ntraj = 500, rng = rng, progress_bar = Val(false))
    @test sol_me_mix.expect ≈ sol_mc_mix.expect atol = 1.0e-2 * length(tlist)

    @testset "Type Inference (mesolve)" begin
        @inferred mesolveProblem(Hq, ψe, tlist_q, (c_op_q,), e_ops = (Pe,), params = pq, progress_bar = Val(false))
        @inferred mesolve(Hq, ψe, tlist_q, (c_op_q,), e_ops = (Pe,), params = pq, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops_td, e_ops = e_ops, params = p, progress_bar = Val(false))
    end

    @testset "Type Inference (mcsolve)" begin
        @inferred mcsolveEnsembleProblem(
            H,
            ψ0,
            tlist,
            c_ops_td,
            ntraj = 5,
            e_ops = e_ops,
            params = p,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred mcsolve(
            H,
            ψ0,
            tlist,
            c_ops_td,
            ntraj = 5,
            e_ops = e_ops,
            params = p,
            progress_bar = Val(false),
            rng = rng,
        )
    end
end

@testset "ComposedOperator support in time evolution" begin
    N = 10  # number of basis states
    a = destroy(N)
    H = a' * a

    ψ0 = basis(N, 9)  # initial state

    γ0 = 0.5
    γ1(p, t) = sqrt(p[1] * exp(-t))
    c_ops_1 = [QobjEvo(2a, γ1)] # This generates a ScaledOperator internally
    c_ops_2 = [QobjEvo(a, γ1) + QobjEvo(a, γ1)] # This generates a ComposedOperator internally

    e_ops = [a' * a]

    tlist = range(0, 10, 100)
    params = [γ0]
    sol_1 = mesolve(H, ψ0, tlist, c_ops_1, e_ops = e_ops; progress_bar = Val(false), params = params, saveat = tlist)
    sol_2 = mesolve(H, ψ0, tlist, c_ops_2, e_ops = e_ops; progress_bar = Val(false), params = params, saveat = tlist)

    @test sol_1.expect[1, :] ≈ sol_2.expect[1, :] atol = 1.0e-10
    @test all(x -> isapprox(x[1].data, x[2].data; atol = 1.0e-10), zip(sol_1.states, sol_2.states))
end
