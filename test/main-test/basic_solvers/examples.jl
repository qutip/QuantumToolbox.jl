using Test
using QuantumToolbox

@testset "Example: Two qubit model" begin
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
    sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = [sp1 * sm1, sp2 * sm2], progress_bar = false) # Here we don't put Val(false) because we want to test the support for Bool type
    sol_mc = mcsolve(H, psi0, t_l, c_ops, e_ops = [sp1 * sm1, sp2 * sm2], progress_bar = Val(false))
    @test sum(abs.(sol_mc.expect[1:2, :] .- sol_me.expect[1:2, :])) / length(t_l) < 0.1
    @test expect(sp1 * sm1, sol_me.states[end]) ≈ expect(sigmap() * sigmam(), ptrace(sol_me.states[end], 1))
end

@testset "Example: Qubit driven by two sequential cosine pulses" begin
    # settings of pulses
    A1 = rand() # Rabi amplitude of pulse 1
    A2 = rand() # Rabi amplitude of pulse 2
    t1 = 0 # starting time of pulse 1
    t2 = 100 # starting time of pulse 2
    T = 10 # duration of each pulse
    pulse1(p, t) = ((t1 <= t) && (t <= t1 + T)) ? (0.5 * p.A1 * π / T) * (1 - cos(2π * (t - t1) / T)) : 0
    pulse2(p, t) = ((t2 <= t) && (t <= t2 + T)) ? (0.5 * p.A2 * π / T) * (1 - cos(2π * (t - t2) / T)) : 0

    # time‐dependent Hamiltonian
    σy = sigmay()
    H = QobjEvo(σy, pulse1) + QobjEvo(σy, pulse2)

    # initial state
    ψ0 = basis(2, 0)
    ρ0 = ket2dm(ψ0) # to force using mesolve

    # sesolve & mesolve
    tlist = range(0.0, stop = 140.0, length = 10000) # don't change tlist, to check the second pulse contribute correctly under 'tstops=tlist' setting
    e_ops = [sigmax(), sigmaz()]
    params = (A1 = A1, A2 = A2)
    prob_se = sesolveProblem(H, ψ0, tlist; e_ops = e_ops, params = params, progress_bar = Val(false))
    prob_me = mesolveProblem(H, ρ0, tlist; e_ops = e_ops, params = params, progress_bar = Val(false))
    sol_se = sesolve(prob_se)
    sol_me = mesolve(prob_me)

    # analytic solution
    function θ(t, t_start, A)
        # a dummy comment to avoid julia formatter
        if t < t_start
            return 0
        elseif (t_start <= t) && (t <= t_start + T)
            return (A * π / T) * (t - t_start - (T * sin(2π * (t - t_start) / T) / (2π)))
        else
            return A * π
        end
    end
    θlist = θ.(tlist, t1, A1) + θ.(tlist, t2, A2)
    X_analytic = sin.(θlist)
    Z_analytic = cos.(θlist)

    @test prob_se.prob.kwargs[:tstops] == tlist # tstops should be equal to tlist for time-dependent cases
    @test prob_me.prob.kwargs[:tstops] == tlist # tstops should be equal to tlist for time-dependent cases
    @test all(isapprox.(X_analytic, sol_se.expect[1, :]; atol = 1.0e-6))
    @test all(isapprox.(X_analytic, sol_me.expect[1, :]; atol = 1.0e-6))
    @test all(isapprox.(Z_analytic, sol_se.expect[2, :]; atol = 1.0e-6))
    @test all(isapprox.(Z_analytic, sol_me.expect[2, :]; atol = 1.0e-6))
end
