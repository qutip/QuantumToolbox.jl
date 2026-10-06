#!PTR_MULTITHREAD
using Test
using QuantumToolbox
import SciMLOperators: ScaledOperator
import Statistics: mean
import Random: MersenneTwister

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "mcsolve" begin
    @test Base.Threads.nthreads() > 1 # make sure the test runs with multi-threading

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
    sol_mc3 = mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), keep_runs_results = Val(true))
    sol_mc_states =
        mcsolve(H, ψ0, tlist, c_ops, saveat = saveat, progress_bar = Val(false), keep_runs_results = Val(true))

    # also test function average_states
    # average the states from all trajectories, and then calculate the expectation value
    expect_mc_states_mean = expect.(Ref(e_ops[1]), average_states(sol_mc_states))

    @test prob_mc.prob.f.f isa ScaledOperator
    # the constant matrices are shared between trajectories, while the scalar coefficients are copied
    L_mc = prob_mc.prob.f.f
    L_mc_copy = @inferred QuantumToolbox._copy_for_trajectory(L_mc)
    @test L_mc_copy.L.A === L_mc.L.A
    @test L_mc_copy.λ !== L_mc.λ
    @test !haskey(prob_mc.prob.kwargs, :tstops) # tstops should not exist for time-independent cases
    @test sum(abs, sol_mc.expect .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, average_expect(sol_mc3) .- sol_me.expect) / length(tlist) < 0.1
    @test sum(abs, expect_mc_states_mean .- vec(sol_me.expect[1, saveat_idxs])) / length(tlist) < 0.1
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

    # check that save_end = false works as expected
    # the states at the end of each trajectory are not saved, but expectation values are still saved
    sol_mc_save_end_false =
        mcsolve(H, ψ0, tlist, c_ops, e_ops = e_ops, save_end = false, progress_bar = Val(false), ntraj = 5)
    @test length(sol_mc_save_end_false.states) == length(sol_mc_save_end_false.times_states) == 0
    @test size(sol_mc_save_end_false.expect) == (length(e_ops), length(tlist))

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
            atol = 1.0e-6,
        ),
    )
    @test average_expect(sol_mc_states) === nothing
    @test std_expect(sol_mc_states) === nothing
    @test_throws ArgumentError std_expect(sol_mc)

    @testset "Memory Allocations (mcsolve)" begin
        ntraj = 100
        for keep_runs_results in (Val(false), Val(true))
            n1 = 125
            n2 = 115

            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                e_ops = e_ops,
                ntraj = ntraj,
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n1 * ntraj + 600 # 125 allocations per trajectory + 600 for initialization

            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            ) # Warm-up
            allocs_tot = @allocations mcsolve(
                H,
                ψ0,
                tlist,
                c_ops,
                ntraj = ntraj,
                saveat = [tlist[end]],
                progress_bar = Val(false),
                keep_runs_results = keep_runs_results,
            )
            @test allocs_tot < n2 * ntraj + 300 # 115 allocations per trajectory + 300 for initialization
        end
    end

    @testset "Type Inference (mcsolve)" begin
        a = TESetup.a
        rng = TESetup.rng

        ens_prob_mc = @inferred mcsolveEnsembleProblem(
            H,
            ψ0,
            tlist,
            c_ops,
            ntraj = 5,
            e_ops = e_ops,
            progress_bar = Val(false),
            rng = rng,
        )
        sol_mc_prob = @inferred mcsolve(ens_prob_mc, keep_runs_results = Val(true)) # ntraj is taken from the problem
        @test sol_mc_prob.ntraj == 5
        @test size(sol_mc_prob.states, 1) == 5
        @inferred mcsolve(H, ψ0, tlist, c_ops, ntraj = 5, e_ops = e_ops, progress_bar = Val(false), rng = rng)
        @inferred mcsolve(H, ψ0, tlist, c_ops, ntraj = 5, progress_bar = Val(true), rng = rng) # test progress bar
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

@testset "mcsolve_map" begin
    # Get parameters from TESetup to simplify the code
    N = TESetup.N
    a = TESetup.a
    σz = TESetup.σz
    σm = TESetup.σm
    c_ops = TESetup.c_ops
    e_ops = TESetup.e_ops
    γ = TESetup.γ

    g = 0.01

    ψ_0_e = tensor(fock(N, 0), basis(2, 0))
    ψ_1_g = tensor(fock(N, 1), basis(2, 1))

    ψ0_list = [ψ_0_e, ψ_1_g]
    ωc_list = [1, 1.01, 1.02]
    ωq_list = [0.96, 0.98]

    tlist = range(0, 10 / γ, 100)

    ωc_fun(p, t) = p[1]
    ωq_fun(p, t) = p[2]
    H = QobjEvo(a' * a, ωc_fun) + QobjEvo(σz / 2, ωq_fun) + g * (a' * σm + a * σm')

    sols_me = mesolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    # Test with multiple initial states but no params
    sols0 = mcsolve_map(TESetup.H, ψ0_list, tlist, c_ops; e_ops = e_ops, ntraj = 10, progress_bar = Val(false))
    # Test with single initial state
    sols1 = mcsolve_map(H, ψ_0_e, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), ntraj = 10, progress_bar = Val(false))
    # Test with multiple initial states
    sols2 = mcsolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    # Test with the states saved for each trajectory
    sols3 = mcsolve_map(
        H,
        ψ0_list,
        tlist,
        c_ops;
        saveat = tlist[(end - 4):end],
        params = (ωc_list, ωq_list),
        ntraj = 10,
        progress_bar = Val(false),
        keep_runs_results = Val(true),
    )

    @test size(sols0) == (2,)
    @test sols0 isa Vector{<:TimeEvolutionMCSol}
    @test size(sols1) == (1, 3, 2)
    @test sols1 isa Array{<:TimeEvolutionMCSol}
    @test size(sols2) == (2, 3, 2)
    @test sols2 isa Array{<:TimeEvolutionMCSol}
    @test size(sols3) == (2, 3, 2)
    @test all(sol -> sol.ntraj == 10, sols1)
    @test all(sol -> sol.ntraj == 500, sols2) # ntraj = 500 by default
    @test all(sol -> size(sol.expect) == (length(e_ops), length(tlist)), sols2)
    @test all(sol -> length(sol.col_times) == length(sol.col_which) == sol.ntraj, sols2)
    @test all(sol -> size(sol.states) == (10, 5), sols3)
    @test all(sol -> sol.expect === nothing, sols3)

    # Compare each combination of initial state and parameters with mesolve_map
    for I in eachindex(sols2)
        @test sum(abs, sols2[I].expect .- sols_me[I].expect) / length(tlist) < 0.1
    end

    # Test with parameter-dependent collapse operators
    γ_fun(p, t) = sqrt(p[1])
    γ_list = [0.05, 0.2]
    sols_γ = mcsolve_map(
        TESetup.H,
        ψ_1_g,
        tlist,
        (QobjEvo(a, γ_fun), sqrt(γ) * σm);
        e_ops = e_ops,
        params = (γ_list,),
        progress_bar = Val(false),
    )
    for (i, γ_a) in enumerate(γ_list)
        sol_me = mesolve(TESetup.H, ψ_1_g, tlist, (sqrt(γ_a) * a, sqrt(γ) * σm); e_ops = e_ops, progress_bar = Val(false))
        @test sum(abs, sols_γ[1, i].expect .- sol_me.expect) / length(tlist) < 0.1
    end

    # Test reproducibility
    sols_rng1 = mcsolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), ntraj = 10, rng = MersenneTwister(1234), progress_bar = Val(false), keep_runs_results = Val(true))
    sols_rng2 = mcsolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), ntraj = 10, rng = MersenneTwister(1234), progress_bar = Val(false), keep_runs_results = Val(true))
    @test all(I -> sols_rng1[I].expect == sols_rng2[I].expect, eachindex(sols_rng1))
    @test all(I -> sols_rng1[I].col_times == sols_rng2[I].col_times, eachindex(sols_rng1))

    @test_throws ArgumentError mcsolve_map(H, ψ0_list, tlist; params = (ωc_list, ωq_list), progress_bar = Val(false))

    @testset "Type Inference mcsolve_map" begin
        @inferred mcsolve_map(TESetup.H, ψ0_list, tlist, c_ops; e_ops = e_ops, ntraj = 5, progress_bar = Val(true)) # no params, but test progress bar
        @inferred mcsolve_map(
            H,
            ψ0_list,
            tlist,
            c_ops;
            e_ops = e_ops,
            params = (ωc_list, ωq_list),
            ntraj = 5,
            progress_bar = Val(false),
        )
        @inferred mcsolve_map(H, ψ0_list, tlist, c_ops; params = (ωc_list, ωq_list), ntraj = 5, progress_bar = Val(false), keep_runs_results = Val(true))
    end
end
