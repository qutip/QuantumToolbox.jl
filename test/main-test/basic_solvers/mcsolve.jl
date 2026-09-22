using Test
using QuantumToolbox
import SciMLOperators: ScaledOperator
import Statistics: mean

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "mcsolve" begin
    

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

    @test prob_mc.prob.f.f isa ScaledOperator
    @test !haskey(prob_mc.prob.kwargs, :tstops) # tstops should not exist for time-independent cases
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
            n1 = 145
            n2 = 135

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
            @test allocs_tot < n1 * ntraj + 600 # 150 allocations per trajectory + 600 for initialization

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