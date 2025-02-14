@testset "Time Evolution and Partial Trace" verbose = true begin
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

    ψ0_int = Qobj(round.(Int, real.(ψ0.data)), dims = ψ0.dims) # Used for testing the type inference

    @testset "sesolve" begin
        tlist = range(0, 20 * 2π / g, 1000)
        saveat_idxs = 500:900
        saveat = tlist[saveat_idxs]

        prob = sesolveProblem(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
        sol = sesolve(prob)
        sol2 = sesolve(H, ψ0, tlist, progress_bar = Val(false))
        sol3 = sesolve(H, ψ0, tlist, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))
        sol_string = sprint((t, s) -> show(t, "text/plain", s), sol)
        sol_string2 = sprint((t, s) -> show(t, "text/plain", s), sol2)

        ## Analytical solution for the expectation value of a' * a
        Ω_rabi = sqrt(g^2 + ((ωc - ωq) / 2)^2)
        amp_rabi = g^2 / Ω_rabi^2
        ##

        @test prob.prob.f.f isa MatrixOperator
        @test sum(abs.(sol.expect[1, :] .- amp_rabi .* sin.(Ω_rabi * tlist) .^ 2)) / length(tlist) < 0.1
        @test length(sol.times) == length(tlist)
        @test length(sol.states) == 1
        @test size(sol.expect) == (length(e_ops), length(tlist))
        @test length(sol2.times) == length(tlist)
        @test length(sol2.states) == length(tlist)
        @test sol2.expect === nothing
        @test length(sol3.times) == length(tlist)
        @test length(sol3.states) == length(saveat)
        @test size(sol3.expect) == (length(e_ops), length(tlist))
        @test sol.expect[1, saveat_idxs] ≈ expect(e_ops[1], sol3.states) atol = 1e-6
        @test sol_string ==
              "Solution of time evolution\n" *
              "(return code: $(sol.retcode))\n" *
              "--------------------------\n" *
              "num_states = $(length(sol.states))\n" *
              "num_expect = $(size(sol.expect, 1))\n" *
              "ODE alg.: $(sol.alg)\n" *
              "abstol = $(sol.abstol)\n" *
              "reltol = $(sol.reltol)\n"
        @test sol_string2 ==
              "Solution of time evolution\n" *
              "(return code: $(sol2.retcode))\n" *
              "--------------------------\n" *
              "num_states = $(length(sol2.states))\n" *
              "num_expect = 0\n" *
              "ODE alg.: $(sol2.alg)\n" *
              "abstol = $(sol2.abstol)\n" *
              "reltol = $(sol2.reltol)\n"

        tlist1 = Float64[]
        tlist2 = [0, 0.2, 0.1]
        tlist3 = [0, 0.1, 0.1, 0.2]
        @test_throws ArgumentError sesolve(H, ψ0, tlist1, progress_bar = Val(false))
        @test_throws ArgumentError sesolve(H, ψ0, tlist2, progress_bar = Val(false))
        @test_throws ArgumentError sesolve(H, ψ0, tlist3, progress_bar = Val(false))

        @testset "Memory Allocations" begin
            allocs_tot = @allocations sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false)) # Warm-up
            allocs_tot = @allocations sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
            @test allocs_tot < 150

            allocs_tot = @allocations sesolve(H, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false)) # Warm-up
            allocs_tot = @allocations sesolve(H, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false))
            @test allocs_tot < 100
        end

        @testset "Type Inference sesolve" begin
            @inferred sesolveProblem(H, ψ0, tlist, progress_bar = Val(false))
            @inferred sesolveProblem(H, ψ0, [0, 10], progress_bar = Val(false))
            @inferred sesolveProblem(H, ψ0_int, tlist, progress_bar = Val(false))
            @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
            @inferred sesolve(H, ψ0, tlist, progress_bar = Val(false))
            @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))
            @inferred sesolve(H, ψ0, tlist, e_ops = (a' * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple of different types
        end
    end

    @testset "mesolve, mcsolve, ssesolve and smesolve" begin
        tlist = range(0, 10 / γ, 100)
        saveat_idxs = 50:90
        saveat = tlist[saveat_idxs]

        prob_me = mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        sol_me = mesolve(prob_me)
        sol_me2 = mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
        sol_me3 = mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))
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
        sol_mc_states = mcsolve(H, ψ0, tlist, c_ops, saveat = saveat, progress_bar = Val(false))
        sol_mc_states2 = mcsolve(
            H,
            ψ0,
            tlist,
            c_ops,
            saveat = saveat,
            progress_bar = Val(false),
            jump_callback = DiscreteLindbladJumpCallback(),
        )
        sol_sse = ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        sol_sme = smesolve(H, ψ0, tlist, c_ops_sme, sc_ops_sme, e_ops = e_ops, progress_bar = Val(false))

        ρt_mc = [ket2dm.(normalize.(states)) for states in sol_mc_states.states]
        expect_mc_states = mapreduce(states -> expect.(Ref(e_ops[1]), states), hcat, ρt_mc)
        expect_mc_states_mean = sum(expect_mc_states, dims = 2) / size(expect_mc_states, 2)

        ρt_mc2 = [ket2dm.(normalize.(states)) for states in sol_mc_states2.states]
        expect_mc_states2 = mapreduce(states -> expect.(Ref(e_ops[1]), states), hcat, ρt_mc2)
        expect_mc_states_mean2 = sum(expect_mc_states2, dims = 2) / size(expect_mc_states2, 2)

        sol_me_string = sprint((t, s) -> show(t, "text/plain", s), sol_me)
        sol_mc_string = sprint((t, s) -> show(t, "text/plain", s), sol_mc)
        sol_mc_string_states = sprint((t, s) -> show(t, "text/plain", s), sol_mc_states)
        sol_sse_string = sprint((t, s) -> show(t, "text/plain", s), sol_sse)
        sol_sme_string = sprint((t, s) -> show(t, "text/plain", s), sol_sme)
        @test prob_me.prob.f.f isa MatrixOperator
        @test prob_mc.prob.f.f isa MatrixOperator
        @test sum(abs, sol_mc.expect .- sol_me.expect) / length(tlist) < 0.1
        @test sum(abs, sol_mc2.expect .- sol_me.expect) / length(tlist) < 0.1
        @test sum(abs, vec(expect_mc_states_mean) .- vec(sol_me.expect[1, saveat_idxs])) / length(tlist) < 0.1
        @test sum(abs, vec(expect_mc_states_mean2) .- vec(sol_me.expect[1, saveat_idxs])) / length(tlist) < 0.1
        @test sum(abs, sol_sse.expect .- sol_me.expect) / length(tlist) < 0.1
        @test sum(abs, sol_sme.expect .- sol_me.expect) / length(tlist) < 0.1
        @test length(sol_me.times) == length(tlist)
        @test length(sol_me.states) == 1
        @test size(sol_me.expect) == (length(e_ops), length(tlist))
        @test length(sol_me2.times) == length(tlist)
        @test length(sol_me2.states) == length(tlist)
        @test sol_me2.expect === nothing
        @test length(sol_me3.times) == length(tlist)
        @test length(sol_me3.states) == length(saveat)
        @test size(sol_me3.expect) == (length(e_ops), length(tlist))
        @test sol_me3.expect[1, saveat_idxs] ≈ expect(e_ops[1], sol_me3.states) atol = 1e-6
        @test length(sol_mc.times) == length(tlist)
        @test size(sol_mc.expect) == (length(e_ops), length(tlist))
        @test length(sol_mc_states.times) == length(tlist)
        @test sol_mc_states.expect === nothing
        @test length(sol_sse.times) == length(tlist)
        @test size(sol_sse.expect) == (length(e_ops), length(tlist))
        @test length(sol_sme.times) == length(tlist)
        @test size(sol_sme.expect) == (length(e_ops), length(tlist))
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
              "num_trajectories = $(sol_mc.ntraj)\n" *
              "num_states = $(length(sol_mc.states[1]))\n" *
              "num_expect = $(size(sol_mc.expect, 1))\n" *
              "ODE alg.: $(sol_mc.alg)\n" *
              "abstol = $(sol_mc.abstol)\n" *
              "reltol = $(sol_mc.reltol)\n"
        @test sol_mc_string_states ==
              "Solution of quantum trajectories\n" *
              "(converged: $(sol_mc_states.converged))\n" *
              "--------------------------------\n" *
              "num_trajectories = $(sol_mc_states.ntraj)\n" *
              "num_states = $(length(sol_mc_states.states[1]))\n" *
              "num_expect = 0\n" *
              "ODE alg.: $(sol_mc_states.alg)\n" *
              "abstol = $(sol_mc_states.abstol)\n" *
              "reltol = $(sol_mc_states.reltol)\n"
        @test sol_sse_string ==
              "Solution of stochastic quantum trajectories\n" *
              "(converged: $(sol_sse.converged))\n" *
              "--------------------------------\n" *
              "num_trajectories = $(sol_sse.ntraj)\n" *
              "num_states = $(length(sol_sse.states[1]))\n" *
              "num_expect = $(size(sol_sse.expect, 1))\n" *
              "SDE alg.: $(sol_sse.alg)\n" *
              "abstol = $(sol_sse.abstol)\n" *
              "reltol = $(sol_sse.reltol)\n"
        @test sol_sme_string ==
              "Solution of stochastic quantum trajectories\n" *
              "(converged: $(sol_sme.converged))\n" *
              "--------------------------------\n" *
              "num_trajectories = $(sol_sme.ntraj)\n" *
              "num_states = $(length(sol_sme.states[1]))\n" *
              "num_expect = $(size(sol_sme.expect, 1))\n" *
              "SDE alg.: $(sol_sme.alg)\n" *
              "abstol = $(sol_sme.abstol)\n" *
              "reltol = $(sol_sme.reltol)\n"

        tlist1 = Float64[]
        tlist2 = [0, 0.2, 0.1]
        tlist3 = [0, 0.1, 0.1, 0.2]
        @test_throws ArgumentError mesolve(H, ψ0, tlist1, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError mesolve(H, ψ0, tlist2, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError mesolve(H, ψ0, tlist3, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError mcsolve(H, ψ0, tlist1, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError mcsolve(H, ψ0, tlist2, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError mcsolve(H, ψ0, tlist3, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError ssesolve(H, ψ0, tlist1, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError ssesolve(H, ψ0, tlist2, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError ssesolve(H, ψ0, tlist3, c_ops, progress_bar = Val(false))
        @test_throws ArgumentError smesolve(H, ψ0, tlist1, c_ops_sme, sc_ops_sme, progress_bar = Val(false))
        @test_throws ArgumentError smesolve(H, ψ0, tlist2, c_ops_sme, sc_ops_sme, progress_bar = Val(false))
        @test_throws ArgumentError smesolve(H, ψ0, tlist3, c_ops_sme, sc_ops_sme, progress_bar = Val(false))

        # Time-Dependent Hamiltonian
        # ssesolve is slow to be run on CI. It is not removed from the test because it may be useful for testing in more powerful machines.

        ωd = 1.02
        F = 0.05

        # Time Evolution in the drive frame

        H_dr_fr = H - ωd * a' * a - ωd * σz / 2 + F * (a + a')

        rng = MersenneTwister(12)

        tlist = range(0, 10 / γ, 1000)

        sol_se = sesolve(H_dr_fr, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
        sol_me = mesolve(H_dr_fr, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        sol_mc = mcsolve(H_dr_fr, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)
        # sol_sse = ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)

        # Time Evolution in the lab frame

        coef1(p, t) = p.F * exp(1im * p.ωd * t)
        coef2(p, t) = p.F * exp(-1im * p.ωd * t)

        H_td = (H, (a, coef1), (a', coef2))
        p = (F = F, ωd = ωd)

        sol_se_td = sesolve(H_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
        sol_me_td = mesolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        sol_mc_td = mcsolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)
        # sol_sse_td = ssesolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)

        @test sol_se.expect ≈ sol_se_td.expect atol = 1e-6 * length(tlist)
        @test sol_me.expect ≈ sol_me_td.expect atol = 1e-6 * length(tlist)
        @test sol_mc.expect ≈ sol_mc_td.expect atol = 1e-2 * length(tlist)
        # @test sol_sse.expect ≈ sol_sse_td.expect atol = 1e-2 * length(tlist)

        H_td2 = QobjEvo(H_td)
        L_td = liouvillian(H_td2)

        sol_se_td2 = sesolve(H_td2, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
        sol_me_td2 = mesolve(L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        sol_mc_td2 = mcsolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)
        # sol_sse_td2 =
        # ssesolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p, rng = rng)

        @test sol_se.expect ≈ sol_se_td2.expect atol = 1e-6 * length(tlist)
        @test sol_me.expect ≈ sol_me_td2.expect atol = 1e-6 * length(tlist)
        @test sol_mc.expect ≈ sol_mc_td2.expect atol = 1e-2 * length(tlist)
        # @test sol_sse.expect ≈ sol_sse_td2.expect atol = 1e-2 * length(tlist)

        @testset "Memory Allocations (mesolve)" begin
            # We predefine the Liouvillian to avoid to count the allocations of the liouvillian function
            L = liouvillian(H, c_ops)
            L_td = QobjEvo((liouvillian(H, c_ops), (liouvillian(a), coef1), (liouvillian(a'), coef2)))

            allocs_tot = @allocations mesolve(L, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false)) # Warm-up
            allocs_tot = @allocations mesolve(L, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false))
            @test allocs_tot < 210

            allocs_tot = @allocations mesolve(L, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false)) # Warm-up
            allocs_tot = @allocations mesolve(L, ψ0, tlist, saveat = [tlist[end]], progress_bar = Val(false))
            @test allocs_tot < 120

            allocs_tot = @allocations mesolve(L_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p) # Warm-up
            allocs_tot = @allocations mesolve(L_td, ψ0, tlist, e_ops = e_ops, progress_bar = Val(false), params = p)
            @test allocs_tot < 210

            allocs_tot =
                @allocations mesolve(L_td, ψ0, tlist, progress_bar = Val(false), saveat = [tlist[end]], params = p) # Warm-up
            allocs_tot =
                @allocations mesolve(L_td, ψ0, tlist, progress_bar = Val(false), saveat = [tlist[end]], params = p)
            @test allocs_tot < 120
        end

        @testset "Memory Allocations (mcsolve)" begin
            ntraj = 100
            allocs_tot =
                @allocations mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, ntraj = ntraj, progress_bar = Val(false)) # Warm-up
            allocs_tot =
                @allocations mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, ntraj = ntraj, progress_bar = Val(false))
            @test allocs_tot < 160 * ntraj + 500 # 150 allocations per trajectory + 500 for initialization

            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            ) # Warm-up
            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            )
            @test allocs_tot < 160 * ntraj + 300 # 100 allocations per trajectory + 300 for initialization
        end

        @testset "Memory Allocations (ssesolve)" begin
            allocs_tot =
                @allocations ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, ntraj = 100, progress_bar = Val(false)) # Warm-up
            allocs_tot =
                @allocations ssesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, ntraj = 100, progress_bar = Val(false))
            @test allocs_tot < 1950000 # TODO: Fix this high number of allocations

            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = 100,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            ) # Warm-up
            allocs_tot = @allocations ssesolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = 100,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            )
            @test allocs_tot < 570000 # TODO: Fix this high number of allocations
        end

        @testset "Memory Allocations (smesolve)" begin
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                e_ops = e_ops,
                ntraj = 100,
                progress_bar = Val(false),
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                e_ops = e_ops,
                ntraj = 100,
                progress_bar = Val(false),
            )
            @test allocs_tot < 2750000 # TODO: Fix this high number of allocations

            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                ntraj = 100,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            ) # Warm-up
            allocs_tot = @allocations smesolve(
                H,
                ψ0,
                tlist,
                c_ops_sme,
                sc_ops_sme,
                ntraj = 100,
                saveat = [tlist[end]],
                progress_bar = Val(false),
            )
            @test allocs_tot < 570000 # TODO: Fix this high number of allocations
        end

        @testset "Type Inference mesolve" begin
            coef(p, t) = exp(-t)
            ad_t = QobjEvo(a', coef)
            @inferred mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolveProblem(H, ψ0, [0, 10], c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolveProblem(H, ψ0_int, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
            @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
            @inferred mesolve(H, ψ0, tlist, (a, ad_t), e_ops = (a' * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple
            @inferred mesolve(H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
            @inferred mesolve(H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
            @inferred mesolve(L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        end

        @testset "Type Inference mcsolve" begin
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
            @inferred mcsolve(H, ψ0_int, tlist, c_ops, ntraj = 5, progress_bar = Val(false), rng = rng)
            @inferred mcsolve(
                H,
                ψ0,
                tlist,
                (a, a'),
                e_ops = (a' * a, a'),
                ntraj = 5,
                progress_bar = Val(false),
                rng = rng,
            ) # We test the type inference for Tuple of different types
            @inferred mcsolve(
                H_td,
                ψ0,
                tlist,
                c_ops,
                ntraj = 5,
                e_ops = e_ops,
                progress_bar = Val(false),
                params = p,
                rng = rng,
            )
        end

        @testset "Type Inference ssesolve" begin
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
            @inferred ssesolve(
                H,
                ψ0,
                tlist,
                c_ops_tuple,
                ntraj = 5,
                e_ops = e_ops,
                progress_bar = Val(false),
                rng = rng,
            )
            @inferred ssesolve(H, ψ0, tlist, c_ops_tuple, ntraj = 5, progress_bar = Val(true), rng = rng)
            @inferred ssesolve(H, ψ0, [0, 10], c_ops_tuple, ntraj = 5, progress_bar = Val(false), rng = rng)
            @inferred ssesolve(H, ψ0_int, tlist, c_ops_tuple, ntraj = 5, progress_bar = Val(false), rng = rng)
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
                H_td,
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

        @testset "Type Inference smesolve" begin
            c_ops_sme_tuple = Tuple(c_ops_sme) # To avoid type instability, we must have a Tuple instead of a Vector
            sc_ops_sme_tuple = Tuple(sc_ops_sme) # To avoid type instability, we must have a Tuple instead of a Vector
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
                ψ0_int,
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

        @testset "mcsolve, ssesolve and smesolve reproducibility" begin
            rng = MersenneTwister(1234)
            sol_mc1 = mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)
            rng = MersenneTwister(1234)
            sol_sse1 = ssesolve(H, ψ0, tlist, c_ops, ntraj = 50, e_ops = e_ops, progress_bar = Val(false), rng = rng)
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
            )

            rng = MersenneTwister(1234)
            sol_mc2 = mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), rng = rng)
            rng = MersenneTwister(1234)
            sol_sse2 = ssesolve(H, ψ0, tlist, c_ops, ntraj = 50, e_ops = e_ops, progress_bar = Val(false), rng = rng)
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
            )

            rng = MersenneTwister(1234)
            sol_mc3 = mcsolve(H, ψ0, tlist, c_ops, ntraj = 510, e_ops = e_ops, progress_bar = Val(false), rng = rng)
            rng = MersenneTwister(1234)
            sol_sse3 = ssesolve(H, ψ0, tlist, c_ops, ntraj = 60, e_ops = e_ops, progress_bar = Val(false), rng = rng)
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
            )

            @test sol_mc1.expect ≈ sol_mc2.expect atol = 1e-10
            @test sol_mc1.runs_expect ≈ sol_mc2.runs_expect atol = 1e-10
            @test sol_mc1.col_times ≈ sol_mc2.col_times atol = 1e-10
            @test sol_mc1.col_which ≈ sol_mc2.col_which atol = 1e-10

            @test sol_mc1.runs_expect ≈ sol_mc3.runs_expect[:, 1:500, :] atol = 1e-10

            @test sol_sse1.expect ≈ sol_sse2.expect atol = 1e-10
            @test sol_sse1.runs_expect ≈ sol_sse2.runs_expect atol = 1e-10

            @test sol_sse1.runs_expect ≈ sol_sse3.runs_expect[:, 1:50, :] atol = 1e-10

            @test sol_sme1.expect ≈ sol_sme2.expect atol = 1e-10
            @test sol_sme1.runs_expect ≈ sol_sme2.runs_expect atol = 1e-10

            @test sol_sme1.runs_expect ≈ sol_sme3.runs_expect[:, 1:50, :] atol = 1e-10
        end
    end

    @testset "exceptions" begin
        N = 10
        a = destroy(N)
        H = a' * a
        c_ops = [sqrt(0.1) * a]
        psi0 = basis(N, 3)
        t_l = LinRange(0, 100, 1000)
        psi_wrong = basis(N - 1, 3)
        @test_throws DimensionMismatch sesolve(H, psi_wrong, t_l)
        @test_throws DimensionMismatch mesolve(H, psi_wrong, t_l, c_ops)
        @test_throws DimensionMismatch mcsolve(H, psi_wrong, t_l, c_ops)
        @test_throws ArgumentError sesolve(H, psi0, t_l, save_idxs = [1, 2])
        @test_throws ArgumentError mesolve(H, psi0, t_l, c_ops, save_idxs = [1, 2])
        @test_throws ArgumentError mcsolve(H, psi0, t_l, c_ops, save_idxs = [1, 2])
    end

    @testset "example" begin
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
end
