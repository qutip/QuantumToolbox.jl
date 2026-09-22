using Test
using QuantumToolbox

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "ssesolve" begin
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
        @inferred ssesolve(H, ψ0, tlist, c_ops_tuple, ntraj = 5, progress_bar = Val(true), rng = rng) # test progress bar
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