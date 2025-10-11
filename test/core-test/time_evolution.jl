@testmodule TESetup begin
    using QuantumToolbox
    using Random

    # Global definition of the system
    N = 10
    a = kron(destroy(N), qeye(2))
    σm = kron(qeye(N), sigmam())
    σz = qeye(N) ⊗ sigmaz()

    g = 0.01
    ωc = 1
    ωq = 0.99
    γ = 0.1
    nth = 0.001

    # Jaynes-Cummings Hamiltonian
    H = ωc * a' * a + ωq / 2 * σz + g * (a' * σm + a * σm')
    ψ0 = kron(fock(N, 0), fock(2, 0))

    e_ops = [a' * a, σz]
    c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a', sqrt(γ * (1 + nth)) * σm, sqrt(γ * nth) * σm']

    sme_η = 0.7 # Efficiency of the homodyne detector for smesolve
    c_ops_sme = [sqrt(1 - sme_η) * op for op in c_ops]
    sc_ops_sme = [sqrt(sme_η) * op for op in c_ops]

    # The following definition is to test the case of `sc_ops` as an `AbstractQuantumObject`
    c_ops_sme2 = c_ops[2:end]
    sc_ops_sme2 = c_ops[1]

    ψ0_int = Qobj(round.(Int, real.(ψ0.data)), dims = ψ0.dims) # Used for testing the type inference

    ψ_wrong = kron(fock(N - 1, 0), fock(2, 0))

    rng = MersenneTwister(12)

    # QobjEvo
    ωd = 1.02
    F = 0.05
    coef1(p, t) = p.F * exp(1im * p.ωd * t)
    coef2(p, t) = p.F * exp(-1im * p.ωd * t)
    p = (F = F, ωd = ωd)
    H_td = (H, (a, coef1), (a', coef2))
    H_td2 = QobjEvo(H_td)
    L_td = liouvillian(H_td2)

    # time list and saveat
    tlist = range(0, 10 / γ, 100)
    saveat_idxs = 50:90
    saveat = tlist[saveat_idxs]

    # time list for testing exceptions
    tlist1 = Float64[]
    tlist2 = [0, 0.2, 0.1]
    tlist3 = [0, 0.1, 0.1, 0.2]

    # mesolve solution used for comparing results from mcsolve, ssesolve, and smesolve with mesolve
    prob_me = mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
    sol_me = mesolve(prob_me)
end

@testitem "sesolve" setup=[TESetup] begin
    using SciMLOperators

    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    e_ops = TESetup.e_ops

    tlist = range(0, 20 * 2π / TESetup.g, 1000)
    saveat_idxs = 500:900
    saveat = tlist[saveat_idxs]

    prob = sesolveProblem(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
    sol = sesolve(prob)
    sol2 = sesolve(H, ψ0, tlist, progress_bar = Val(false))
    sol3 = sesolve(H, ψ0, tlist, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))

    ## Analytical solution for the expectation value of a' * a
    Ω_rabi = sqrt(TESetup.g^2 + ((TESetup.ωc - TESetup.ωq) / 2)^2)
    amp_rabi = TESetup.g^2 / Ω_rabi^2
    ##

    @test prob.prob.f.f isa MatrixOperator
    @test sum(abs.(sol.expect[1, :] .- amp_rabi .* sin.(Ω_rabi * tlist) .^ 2)) / length(tlist) < 0.1
    @test length(sol.times) == length(tlist)
    @test length(sol.times_states) == 1
    @test length(sol.states) == 1
    @test size(sol.expect) == (length(e_ops), length(tlist))
    @test length(sol2.times) == length(tlist)
    @test length(sol2.times_states) == length(tlist)
    @test length(sol2.states) == length(tlist)
    @test sol2.expect === nothing
    @test length(sol3.times) == length(tlist)
    @test length(sol3.times_states) == length(saveat)
    @test length(sol3.states) == length(saveat)
    @test size(sol3.expect) == (length(e_ops), length(tlist))
    @test sol.expect[1, saveat_idxs] ≈ expect(e_ops[1], sol3.states) atol = 1e-6

    sol_string = sprint((t, s) -> show(t, "text/plain", s), sol)
    @test sol_string ==
          "Solution of time evolution\n" *
          "(return code: $(sol.retcode))\n" *
          "--------------------------\n" *
          "num_states = $(length(sol.states))\n" *
          "num_expect = $(size(sol.expect, 1))\n" *
          "ODE alg.: $(sol.alg)\n" *
          "abstol = $(sol.abstol)\n" *
          "reltol = $(sol.reltol)\n"

    sol_string2 = sprint((t, s) -> show(t, "text/plain", s), sol2)
    @test sol_string2 ==
          "Solution of time evolution\n" *
          "(return code: $(sol2.retcode))\n" *
          "--------------------------\n" *
          "num_states = $(length(sol2.states))\n" *
          "num_expect = 0\n" *
          "ODE alg.: $(sol2.alg)\n" *
          "abstol = $(sol2.abstol)\n" *
          "reltol = $(sol2.reltol)\n"

    @test_throws ArgumentError sesolve(H, ψ0, TESetup.tlist1, progress_bar = Val(false))
    @test_throws ArgumentError sesolve(H, ψ0, TESetup.tlist2, progress_bar = Val(false))
    @test_throws ArgumentError sesolve(H, ψ0, TESetup.tlist3, progress_bar = Val(false))
    @test_throws ArgumentError sesolve(H, ψ0, tlist, save_idxs = [1, 2], progress_bar = Val(false))
    @test_throws DimensionMismatch sesolve(H, TESetup.ψ_wrong, tlist, progress_bar = Val(false))

    @testset "Memory Allocations" begin
        allocs_tot = @allocations sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false)) # Warm-up
        allocs_tot = @allocations sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
        @test allocs_tot < 110

        allocs_tot = @allocations sesolve(H, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false)) # Warm-up
        allocs_tot = @allocations sesolve(H, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false))
        @test allocs_tot < 90
    end

    @testset "Type Inference sesolve" begin
        @inferred sesolveProblem(H, ψ0, tlist, progress_bar = Val(false))
        @inferred sesolveProblem(H, ψ0, [0, 10], progress_bar = Val(false))
        @inferred sesolveProblem(H, TESetup.ψ0_int, tlist, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = (TESetup.a' * TESetup.a, TESetup.a'), progress_bar = Val(false)) # We test the type inference for Tuple of different types
    end
end

@testitem "sesolve_map" setup=[TESetup] begin

    # Get parameters from TESetup to simplify the code
    N = TESetup.N
    a = TESetup.a
    σz = TESetup.σz
    σm = TESetup.σm
    e_ops = TESetup.e_ops

    g = 0.01

    ψ_0_e = tensor(fock(N, 0), basis(2, 0))
    ψ_1_g = tensor(fock(N, 1), basis(2, 1))

    ψ0_list = [ψ_0_e, ψ_1_g]
    ωc_list = [1, 1.01, 1.02]
    ωq_list = [0.96, 0.97, 0.98, 0.99]

    tlist = range(0, 20 * 2π / g, 1000)

    ωc_fun(p, t) = p[1]
    ωq_fun(p, t) = p[2]
    H = QobjEvo(a' * a, ωc_fun) + QobjEvo(σz / 2, ωq_fun) + g * (a' * σm + a * σm')

    sols1 = sesolve_map(H, ψ_0_e, tlist; e_ops = e_ops, params = (ωc_list, ωq_list))
    sols2 = sesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    @test size(sols1) == (1, 3, 4)
    @test sols1 isa Array{<:TimeEvolutionSol}
    @test size(sols2) == (2, 3, 4)
    @test sols2 isa Array{<:TimeEvolutionSol}
    for (i, ωc) in enumerate(ωc_list)
        for (j, ωq) in enumerate(ωq_list)
            sol_0_e = sols2[1, i, j]
            sol_1_g = sols2[2, i, j]

            ## Analytical solution for the expectation value of a' * a
            Ω_rabi = sqrt(g^2 + ((ωc - ωq) / 2)^2)
            amp_rabi = g^2 / Ω_rabi^2

            @test sol_0_e.expect[1, :] ≈ amp_rabi .* sin.(Ω_rabi * tlist) .^ 2 atol = 1e-2
            @test sol_1_g.expect[1, :] ≈ 1 .- amp_rabi .* sin.(Ω_rabi * tlist) .^ 2 atol = 1e-2
        end
    end

    @testset "Type Inference sesolve_map" begin
        @inferred sesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    end
end

@testitem "mesolve" setup=[TESetup] begin
    using SciMLOperators

    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    tlist = TESetup.tlist
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    saveat = TESetup.saveat
    sol_me = TESetup.sol_me

    sol_me2 = mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
    sol_me3 = mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))

    # For testing the `OperatorKet` input
    sol_me4 = mesolve(H, operator_to_vector(ket2dm(ψ0)), tlist, c_ops, saveat = saveat, progress_bar = Val(false))

    # Redirect to `sesolve`
    sol_me5 = mesolve(H, ψ0, tlist, progress_bar = Val(false))

    @test TESetup.prob_me.prob.f.f isa MatrixOperator
    @test isket(sol_me5.states[1])
    @test length(sol_me.times) == length(tlist)
    @test length(sol_me.times_states) == 1
    @test length(sol_me.states) == 1
    @test size(sol_me.expect) == (length(e_ops), length(tlist))
    @test length(sol_me2.times) == length(tlist)
    @test length(sol_me2.times_states) == length(tlist)
    @test length(sol_me2.states) == length(tlist)
    @test sol_me2.expect === nothing
    @test length(sol_me3.times) == length(tlist)
    @test length(sol_me3.times_states) == length(saveat)
    @test length(sol_me3.states) == length(saveat)
    @test size(sol_me3.expect) == (length(e_ops), length(tlist))
    @test sol_me3.expect[1, TESetup.saveat_idxs] ≈ expect(e_ops[1], sol_me3.states) atol = 1e-6
    @test all([sol_me3.states[i] ≈ vector_to_operator(sol_me4.states[i]) for i in eachindex(saveat)])

    sol_me_string = sprint((t, s) -> show(t, "text/plain", s), sol_me)
    @test sol_me_string ==
          "Solution of time evolution\n" *
          "(return code: $(sol_me.retcode))\n" *
          "--------------------------\n" *
          "num_states = $(length(sol_me.states))\n" *
          "num_expect = $(size(sol_me.expect, 1))\n" *
          "ODE alg.: $(sol_me.alg)\n" *
          "abstol = $(sol_me.abstol)\n" *
          "reltol = $(sol_me.reltol)\n"

    @test_throws ArgumentError mesolve(H, ψ0, TESetup.tlist1, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mesolve(H, ψ0, TESetup.tlist2, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mesolve(H, ψ0, TESetup.tlist3, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mesolve(H, ψ0, tlist, c_ops, save_idxs = [1, 2], progress_bar = Val(false))
    @test_throws DimensionMismatch mesolve(H, TESetup.ψ_wrong, tlist, c_ops, progress_bar = Val(false))

    @testset "Memory Allocations (mesolve)" begin
        a = TESetup.a
        p = TESetup.p

        # We predefine the Liouvillian to avoid to count the allocations of the liouvillian function
        L = liouvillian(H, c_ops)
        L_td = QobjEvo((liouvillian(H, c_ops), (liouvillian(a), TESetup.coef1), (liouvillian(a'), TESetup.coef2)))

        allocs_tot = @allocations mesolve(L, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false)) # Warm-up
        allocs_tot = @allocations mesolve(L, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
        @test allocs_tot < 180

        allocs_tot = @allocations mesolve(L, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false)) # Warm-up
        allocs_tot = @allocations mesolve(L, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false))
        @test allocs_tot < 110

        allocs_tot = @allocations mesolve(L_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p) # Warm-up
        allocs_tot = @allocations mesolve(L_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
        @test allocs_tot < 180

        allocs_tot = @allocations mesolve(L_td, ψ0, tlist, progress_bar = Val(false), saveat = [tlist[end]], params = p) # Warm-up
        allocs_tot = @allocations mesolve(L_td, ψ0, tlist, progress_bar = Val(false), saveat = [tlist[end]], params = p)
        @test allocs_tot < 110
    end

    @testset "Type Inference (mesolve)" begin
        a = TESetup.a
        p = TESetup.p

        coef(p, t) = exp(-t)
        ad_t = QobjEvo(a', coef)
        @inferred mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolveProblem(H, ψ0, [0, 10], c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolveProblem(H, TESetup.ψ0_int, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, (a, ad_t), e_ops = (a' * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple
        @inferred mesolve(TESetup.H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        @inferred mesolve(TESetup.H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        @inferred mesolve(TESetup.L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
    end
end

@testitem "mesolve_map" setup=[TESetup] begin

    # Get parameters from TESetup to simplify the code
    N = TESetup.N
    a = TESetup.a
    σz = TESetup.σz
    σm = TESetup.σm
    ψ0 = TESetup.ψ0
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    γ = TESetup.γ
    nth = TESetup.nth

    g = 0.01

    ψ_0_e = tensor(fock(N, 0), basis(2, 0))
    ψ_1_g = tensor(fock(N, 1), basis(2, 1))

    ψ0_list = [ψ_0_e, ψ_1_g]
    ωc_list = [1, 1.01, 1.02]
    ωq_list = [0.96, 0.97, 0.98, 0.99]

    tlist = range(0, 10 / γ, 100)

    ωc_fun(p, t) = p[1]
    ωq_fun(p, t) = p[2]
    H = QobjEvo(a' * a, ωc_fun) + QobjEvo(σz / 2, ωq_fun) + g * (a' * σm + a * σm')

    # Test with single initial state
    sols1 = mesolve_map(H, ψ_0_e, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list))
    # Test with multiple initial states
    sols2 = mesolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    # Test redirect to sesolve_map when c_ops is nothing
    sols3 = mesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    @test size(sols1) == (1, 3, 4)
    @test sols1 isa Array{<:TimeEvolutionSol}
    @test size(sols2) == (2, 3, 4)
    @test sols2 isa Array{<:TimeEvolutionSol}
    @test size(sols3) == (2, 3, 4)
    @test sols3 isa Array{<:TimeEvolutionSol}

    # Verify that solutions make physical sense
    for (i, ωc) in enumerate(ωc_list)
        for (j, ωq) in enumerate(ωq_list)
            sol_0_e = sols2[1, i, j]
            sol_1_g = sols2[2, i, j]

            # Check that expectation values are bounded and physical (take real part for physical observables)
            @test all(x -> real(x) >= -1e-4, sol_0_e.expect[1, :]) # a'a should be non-negative (with small tolerance)
            @test all(x -> real(x) >= -1e-4, sol_1_g.expect[1, :])
        end
    end

    # Test with OperatorKet input
    ρ0 = operator_to_vector(ket2dm(ψ_0_e))
    ρ0_list = [operator_to_vector(ket2dm(ψ_0_e)), operator_to_vector(ket2dm(ψ_1_g))]
    sols4 = mesolve_map(H, ρ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    @test size(sols4) == (2, 3, 4)
    @test all(isoperket.(getfield.(sols4, :states) .|> first))

    # Test with Operator input (density matrix)
    dm0_list = [ket2dm(ψ_0_e), ket2dm(ψ_1_g)]
    sols5 =
        mesolve_map(H, dm0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    @test size(sols5) == (2, 3, 4)
    @test sols5 isa Array{<:TimeEvolutionSol}

    @testset "Type Inference mesolve_map" begin
        @inferred mesolve_map(
            H,
            ψ0_list,
            tlist,
            c_ops;
            e_ops = e_ops,
            params = (ωc_list, ωq_list),
            progress_bar = Val(false),
        )
    end
end

@testitem "mcsolve" setup=[TESetup] begin
    using SciMLOperators
    using Statistics

    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    tlist = TESetup.tlist
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    saveat = TESetup.saveat
    saveat_idxs = TESetup.saveat_idxs
    sol_me = TESetup.sol_me

    prob_mc = mcsolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
    sol_mc = mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
    sol_mc2 = mcsolve(
        H,
        ψ0,
        tlist,
        c_ops,
        e_ops = e_ops,
        progress_bar = Val(false),
        jump_callback = DiscreteLindbladJumpCallback(),
    )
    sol_mc3 = mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), keep_runs_results = Val(true))
    sol_mc_states =
        mcsolve(H, ψ0, tlist, c_ops, saveat = saveat, progress_bar = Val(false), keep_runs_results = Val(true))
    sol_mc_states2 = mcsolve(
        H,
        ψ0,
        tlist,
        c_ops,
        saveat = saveat,
        progress_bar = Val(false),
        jump_callback = DiscreteLindbladJumpCallback(),
        keep_runs_results = Val(true),
    )

    # also test function average_states
    # average the states from all trajectories, and then calculate the expectation value
    expect_mc_states_mean = expect.(Ref(e_ops[1]), average_states(sol_mc_states))
    expect_mc_states_mean2 = expect.(Ref(e_ops[1]), average_states(sol_mc_states2))

    @test prob_mc.prob.f.f isa MatrixOperator
    @test sum(abs, sol_mc.expect .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, sol_mc2.expect .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, average_expect(sol_mc3) .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, expect_mc_states_mean .- vec(sol_me.expect[1, saveat_idxs])) / length(tlist) < 0.1
    @test sum(abs, expect_mc_states_mean2 .- vec(sol_me.expect[1, saveat_idxs])) / length(tlist) < 0.1
    @test length(sol_mc.times) == length(tlist)
    @test length(sol_mc.times_states) == 1
    @test size(sol_mc.expect) == (length(e_ops), length(tlist))
    @test size(sol_mc.states) == (1,)
    @test length(sol_mc3.times) == length(tlist)
    @test length(sol_mc3.times_states) == 1
    @test size(sol_mc3.expect) == (length(e_ops), 500, length(tlist)) # ntraj = 500
    @test size(sol_mc3.states) == (500, 1) # ntraj = 500
    @test length(sol_mc_states.times) == length(tlist)
    @test length(sol_mc_states.times_states) == length(saveat)
    @test size(sol_mc_states.states) == (500, length(saveat)) # ntraj = 500
    @test sol_mc_states.expect === nothing

    sol_mc_string = sprint((t, s) -> show(t, "text/plain", s), sol_mc)
    sol_mc_string_states = sprint((t, s) -> show(t, "text/plain", s), sol_mc_states)
    @test sol_mc_string ==
          "Solution of quantum trajectories\n" *
          "(converged: $(sol_mc.converged))\n" *
          "--------------------------------\n" *
          "num_trajectories = $(sol_mc.ntraj)\n" *
          "num_states = $(size(sol_mc.states, ndims(sol_mc.states)))\n" *
          "num_expect = $(size(sol_mc.expect, 1))\n" *
          "ODE alg.: $(sol_mc.alg)\n" *
          "abstol = $(sol_mc.abstol)\n" *
          "reltol = $(sol_mc.reltol)\n"
    @test sol_mc_string_states ==
          "Solution of quantum trajectories\n" *
          "(converged: $(sol_mc_states.converged))\n" *
          "--------------------------------\n" *
          "num_trajectories = $(sol_mc_states.ntraj)\n" *
          "num_states = $(size(sol_mc_states.states, ndims(sol_mc_states.states)))\n" *
          "num_expect = 0\n" *
          "ODE alg.: $(sol_mc_states.alg)\n" *
          "abstol = $(sol_mc_states.abstol)\n" *
          "reltol = $(sol_mc_states.reltol)\n"

    @test_throws ArgumentError mcsolve(H, ψ0, TESetup.tlist1, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mcsolve(H, ψ0, TESetup.tlist2, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mcsolve(H, ψ0, TESetup.tlist3, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError mcsolve(H, ψ0, tlist, c_ops, save_idxs = [1, 2], progress_bar = Val(false))
    @test_throws DimensionMismatch mcsolve(H, TESetup.ψ_wrong, tlist, c_ops, progress_bar = Val(false))

    # test average_states, average_expect, and std_expect
    expvals_all = sol_mc3.expect[:, :, 2:end] # ignore testing initial time point since its standard deviation is a very small value (basically zero)
    stdvals = std_expect(sol_mc3)
    @test average_states(sol_mc) == sol_mc.states
    @test average_expect(sol_mc) == sol_mc.expect
    @test size(stdvals) == (length(e_ops), length(tlist))
    @test all(
        isapprox.(
            stdvals[:, 2:end], # ignore testing initial time point since its standard deviation is a very small value (basically zero)
            dropdims(sqrt.(mean(abs2.(expvals_all), dims = 2) .- abs2.(mean(expvals_all, dims = 2))), dims = 2);
            atol = 1e-6,
        ),
    )
    @test average_expect(sol_mc_states) === nothing
    @test std_expect(sol_mc_states) === nothing
    @test_throws ArgumentError std_expect(sol_mc)

    @testset "Memory Allocations (mcsolve)" begin
        ntraj = 100
        for keep_runs_results in (Val(false), Val(true))
            n1 = QuantumToolbox.getVal(keep_runs_results) ? 120 : 140
            n2 = QuantumToolbox.getVal(keep_runs_results) ? 110 : 130

            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = Val(true),
            ) # Warm-up
            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = Val(true),
            )
            @test allocs_tot < n1 * ntraj + 400 # 150 allocations per trajectory + 500 for initialization

            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = Val(true),
            ) # Warm-up
            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = Val(true),
            )
            @test allocs_tot < n2 * ntraj + 300 # 100 allocations per trajectory + 300 for initialization
        end
    end

    @testset "Type Inference (mcsolve)" begin
        a = TESetup.a
        rng = TESetup.rng

        @inferred mcsolveEnsembleProblem(
            H,
            ψ0,
            tlist,
            c_ops,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred mcsolve(H, ψ0, tlist, c_ops, ntraj = 5, e_ops = e_ops, progress_bar = Val(false), rng = rng)
        @inferred mcsolve(H, ψ0, tlist, c_ops, ntraj = 5, progress_bar = Val(true), rng = rng)
        @inferred mcsolve(H, ψ0, [0, 10], c_ops, ntraj = 5, progress_bar = Val(false), rng = rng)
        @inferred mcsolve(H, TESetup.ψ0_int, tlist, c_ops, ntraj = 5, progress_bar = Val(false), rng = rng)
        @inferred mcsolve(H, ψ0, tlist, (a, a'), e_ops = (a' * a, a'), ntraj = 5, progress_bar = Val(false), rng = rng) # We test the type inference for Tuple of different types
        @inferred mcsolve(
            TESetup.H_td,
            ψ0,
            tlist,
            c_ops,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            params = TESetup.p,
            rng = rng,
        )
    end
end

@testitem "ssesolve" setup=[TESetup] begin
    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    tlist = TESetup.tlist
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    sol_me = TESetup.sol_me

    sol_sse = ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
    sol_sse2 = ssesolve(
        H,
        ψ0,
        tlist,
        c_ops,
        e_ops = e_ops,
        ntraj = 20,
        progress_bar = Val(false),
        store_measurement = Val(true),
    )

    @test sum(abs, sol_sse.expect .- sol_me.expect) / length(tlist) < 0.1
    @test length(sol_sse.times) == length(tlist)
    @test length(sol_sse.times_states) == 1
    @test size(sol_sse.states) == (1,) # ntraj = 500 but keep_runs_results = Val(false)
    @test size(sol_sse.expect) == (length(e_ops), length(tlist))
    @test isnothing(sol_sse.measurement)
    @test size(sol_sse2.measurement) == (length(c_ops), 20, length(tlist) - 1)

    sol_sse_string = sprint((t, s) -> show(t, "text/plain", s), sol_sse)
    @test sol_sse_string ==
          "Solution of stochastic quantum trajectories\n" *
          "(converged: $(sol_sse.converged))\n" *
          "--------------------------------\n" *
          "num_trajectories = $(sol_sse.ntraj)\n" *
          "num_states = $(size(sol_sse.states, ndims(sol_sse.states)))\n" *
          "num_expect = $(size(sol_sse.expect, 1))\n" *
          "SDE alg.: $(sol_sse.alg)\n" *
          "abstol = $(sol_sse.abstol)\n" *
          "reltol = $(sol_sse.reltol)\n"

    @test_throws ArgumentError ssesolve(H, ψ0, TESetup.tlist1, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError ssesolve(H, ψ0, TESetup.tlist2, c_ops, progress_bar = Val(false))
    @test_throws ArgumentError ssesolve(H, ψ0, TESetup.tlist3, c_ops, progress_bar = Val(false))

    @testset "Memory Allocations (ssesolve)" begin
        ntraj = 100
        for keep_runs_results in (Val(false), Val(true))
            n1 = QuantumToolbox.getVal(keep_runs_results) ? 1100 : 1120
            n2 = QuantumToolbox.getVal(keep_runs_results) ? 1000 : 1020

            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n1 * ntraj + 400 # TODO: Fix this high number of allocations

            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n2 * ntraj + 300 # TODO: Fix this high number of allocations
        end
    end

    @testset "Type Inference (ssesolve)" begin
        a = TESetup.a
        rng = TESetup.rng
        p = TESetup.p

        c_ops_tuple = Tuple(c_ops) # To avoid type instability, we must have a Tuple instead of a Vector
        @inferred ssesolveEnsembleProblem(
            H,
            ψ0,
            tlist,
            c_ops_tuple,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred ssesolve(H, ψ0, tlist, c_ops_tuple, ntraj = 5, e_ops = e_ops, progress_bar = Val(false), rng = rng)
        @inferred ssesolve(H, ψ0, tlist, c_ops_tuple, ntraj = 5, progress_bar = Val(true), rng = rng)
        @inferred ssesolve(H, ψ0, [0, 10], c_ops_tuple, ntraj = 5, progress_bar = Val(false), rng = rng)
        @inferred ssesolve(H, TESetup.ψ0_int, tlist, c_ops_tuple, ntraj = 5, progress_bar = Val(false), rng = rng)
        @inferred ssesolve(
            H,
            ψ0,
            tlist,
            c_ops_tuple,
            ntraj = 5,
            e_ops = (a' * a, a'),
            progress_bar = Val(false),
            rng = rng,
        ) # We test the type inference for Tuple of different types
        @inferred ssesolve(
            TESetup.H_td,
            ψ0,
            tlist,
            c_ops_tuple,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            params = p,
            rng = rng,
        )
    end
end

@testitem "smesolve" setup=[TESetup] begin
    using Random

    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    tlist = TESetup.tlist
    c_ops_sme = TESetup.c_ops_sme
    sc_ops_sme = TESetup.sc_ops_sme
    c_ops_sme2 = TESetup.c_ops_sme2
    sc_ops_sme2 = TESetup.sc_ops_sme2
    e_ops = TESetup.e_ops
    sol_me = TESetup.sol_me
    saveat = TESetup.saveat

    sol_sme = smesolve(H, ψ0, tlist, c_ops_sme, sc_ops_sme, e_ops = e_ops, progress_bar = Val(false))
    sol_sme2 = smesolve(
        H,
        ψ0,
        tlist,
        c_ops_sme,
        sc_ops_sme,
        e_ops = e_ops,
        ntraj = 20,
        progress_bar = Val(false),
        store_measurement = Val(true),
    )
    sol_sme3 = smesolve(H, ψ0, tlist, c_ops_sme2, sc_ops_sme2, e_ops = e_ops, progress_bar = Val(false))

    # For testing the `OperatorKet` input
    sol_sme4 = smesolve(
        H,
        ψ0,
        tlist,
        c_ops_sme,
        sc_ops_sme,
        saveat = saveat,
        ntraj = 10,
        progress_bar = Val(false),
        rng = MersenneTwister(12),
    )
    sol_sme5 = smesolve(
        H,
        operator_to_vector(ket2dm(ψ0)),
        tlist,
        c_ops_sme,
        sc_ops_sme,
        saveat = saveat,
        ntraj = 10,
        progress_bar = Val(false),
        rng = MersenneTwister(12),
    )

    @test sum(abs, sol_sme.expect .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, sol_sme3.expect .- sol_me.expect) / length(tlist) < 0.1
    @test length(sol_sme.times) == length(tlist)
    @test length(sol_sme.times_states) == 1
    @test size(sol_sme.states) == (1,) # ntraj = 500 but keep_runs_results = Val(false)
    @test size(sol_sme.expect) == (length(e_ops), length(tlist))
    @test isnothing(sol_sme.measurement)
    @test size(sol_sme2.measurement) == (length(sc_ops_sme), 20, length(tlist) - 1)
    @test all([sol_sme4.states[i] ≈ vector_to_operator(sol_sme5.states[i]) for i in eachindex(saveat)])

    sol_sme_string = sprint((t, s) -> show(t, "text/plain", s), sol_sme)
    @test sol_sme_string ==
          "Solution of stochastic quantum trajectories\n" *
          "(converged: $(sol_sme.converged))\n" *
          "--------------------------------\n" *
          "num_trajectories = $(sol_sme.ntraj)\n" *
          "num_states = $(size(sol_sme.states, ndims(sol_sme.states)))\n" *
          "num_expect = $(size(sol_sme.expect, 1))\n" *
          "SDE alg.: $(sol_sme.alg)\n" *
          "abstol = $(sol_sme.abstol)\n" *
          "reltol = $(sol_sme.reltol)\n"

    @test_throws ArgumentError smesolve(H, ψ0, TESetup.tlist1, c_ops_sme, sc_ops_sme, progress_bar = Val(false))
    @test_throws ArgumentError smesolve(H, ψ0, TESetup.tlist2, c_ops_sme, sc_ops_sme, progress_bar = Val(false))
    @test_throws ArgumentError smesolve(H, ψ0, TESetup.tlist3, c_ops_sme, sc_ops_sme, progress_bar = Val(false))

    @testset "Memory Allocations (smesolve)" begin
        ntraj = 100
        for keep_runs_results in (Val(false), Val(true))
            n1 = QuantumToolbox.getVal(keep_runs_results) ? 1100 : 1120
            n2 = QuantumToolbox.getVal(keep_runs_results) ? 1000 : 1020
            n3 = QuantumToolbox.getVal(keep_runs_results) ? 600 : 620
            n4 = QuantumToolbox.getVal(keep_runs_results) ? 550 : 570

            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n1 * ntraj + 2300 # TODO: Fix this high number of allocations

            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n2 * ntraj + 1500 # TODO: Fix this high number of allocations

            # Diagonal Noise Case
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme2,
                sc_ops_sme2,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme2,
                sc_ops_sme2,
                e_ops = e_ops,
                ntraj = 1,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n3 * ntraj + 1400 # TODO: Fix this high number of allocations

            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme2,
                sc_ops_sme2,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme2,
                sc_ops_sme2,
                ntraj = 1,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n4 * ntraj + 1000 # TODO: Fix this high number of allocations
        end
    end

    @testset "Type Inference (smesolve)" begin
        a = TESetup.a
        rng = TESetup.rng

        # To avoid type instability, we must have a Tuple instead of a Vector
        c_ops_sme_tuple = Tuple(c_ops_sme)
        sc_ops_sme_tuple = Tuple(sc_ops_sme)
        c_ops_sme2_tuple = Tuple(c_ops_sme2)
        sc_ops_sme2_tuple = sc_ops_sme2 # This is an `AbstractQuantumObject`
        @inferred smesolveEnsembleProblem(
            H,
            ψ0,
            tlist,
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred smesolve(
            H,
            ψ0,
            tlist,
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred smesolve(
            H,
            ψ0,
            tlist,
            c_ops_sme2_tuple,
            sc_ops_sme2_tuple,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred smesolve(
            H,
            ψ0,
            tlist,
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            progress_bar = Val(true),
            rng = rng,
        )
        @inferred smesolve(
            H,
            ψ0,
            [0, 10],
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred smesolve(
            H,
            TESetup.ψ0_int,
            tlist,
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            progress_bar = Val(false),
            rng = rng,
        )
        @inferred smesolve(
            H,
            ψ0,
            tlist,
            c_ops_sme_tuple,
            sc_ops_sme_tuple,
            ntraj = 5,
            e_ops = (a' * a, a'),
            progress_bar = Val(false),
            rng = rng,
        ) # We test the type inference for Tuple of different types
    end
end

@testitem "Time-dependent Hamiltonian" setup=[TESetup] begin
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

    @test sol_se.expect ≈ sol_se_td.expect atol = 1e-6 * length(tlist)
    @test sol_me.expect ≈ sol_me_td.expect atol = 1e-6 * length(tlist)
    @test sol_mc.expect ≈ sol_mc_td.expect atol = 1e-2 * length(tlist)
    # @test sol_sse.expect ≈ sol_sse_td.expect atol = 1e-2 * length(tlist)

    sol_se_td2 = sesolve(H_td2, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_me_td2 = mesolve(L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
    sol_mc_td2 = mcsolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)
    # sol_sse_td2 =
    # ssesolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)

    @test sol_se.expect ≈ sol_se_td2.expect atol = 1e-6 * length(tlist)
    @test sol_me.expect ≈ sol_me_td2.expect atol = 1e-6 * length(tlist)
    @test sol_mc.expect ≈ sol_mc_td2.expect atol = 1e-2 * length(tlist)
    # @test sol_sse.expect ≈ sol_sse_td2.expect atol = 1e-2 * length(tlist)
end

@testitem "mcsolve, ssesolve and smesolve reproducibility" setup=[TESetup] begin
    using Random

    # Get parameters from TESetup to simplify the code
    H = TESetup.H
    ψ0 = TESetup.ψ0
    tlist = TESetup.tlist
    c_ops = TESetup.c_ops
    c_ops_sme = TESetup.c_ops_sme
    sc_ops_sme = TESetup.sc_ops_sme
    e_ops = TESetup.e_ops
    rng = TESetup.rng

    rng = MersenneTwister(1234)
    sol_mc1 =
        mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng, keep_runs_results = Val(true))
    rng = MersenneTwister(1234)
    sol_sse1 = ssesolve(
        H,
        ψ0,
        tlist,
        c_ops,
        ntraj = 50,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )
    rng = MersenneTwister(1234)
    sol_sme1 = smesolve(
        H,
        ψ0,
        tlist,
        c_ops_sme,
        sc_ops_sme,
        ntraj = 50,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )

    rng = MersenneTwister(1234)
    sol_mc2 =
        mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng, keep_runs_results = Val(true))
    rng = MersenneTwister(1234)
    sol_sse2 = ssesolve(
        H,
        ψ0,
        tlist,
        c_ops,
        ntraj = 50,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )
    rng = MersenneTwister(1234)
    sol_sme2 = smesolve(
        H,
        ψ0,
        tlist,
        c_ops_sme,
        sc_ops_sme,
        ntraj = 50,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )

    rng = MersenneTwister(1234)
    sol_mc3 = mcsolve(
        H,
        ψ0,
        tlist,
        c_ops,
        ntraj = 510,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )
    rng = MersenneTwister(1234)
    sol_sse3 = ssesolve(
        H,
        ψ0,
        tlist,
        c_ops,
        ntraj = 60,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )
    rng = MersenneTwister(1234)
    sol_sme3 = smesolve(
        H,
        ψ0,
        tlist,
        c_ops_sme,
        sc_ops_sme,
        ntraj = 60,
        e_ops = e_ops,
        progress_bar = Val(false),
        rng = rng,
        keep_runs_results = Val(true),
    )

    @test sol_mc1.expect ≈ sol_mc2.expect atol = 1e-10
    @test sol_mc1.col_times ≈ sol_mc2.col_times atol = 1e-10
    @test sol_mc1.col_which ≈ sol_mc2.col_which atol = 1e-10

    @test sol_mc1.expect ≈ sol_mc3.expect[:, 1:500, :] atol = 1e-10

    @test sol_sse1.expect ≈ sol_sse2.expect atol = 1e-10

    @test sol_sse1.expect ≈ sol_sse3.expect[:, 1:50, :] atol = 1e-10

    @test sol_sme1.expect ≈ sol_sme2.expect atol = 1e-10

    @test sol_sme1.expect ≈ sol_sme3.expect[:, 1:50, :] atol = 1e-10
end

@testitem "example" begin
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
