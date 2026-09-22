using Test
using QuantumToolbox
import Random: MersenneTwister

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "smesolve" begin

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
        ) # test progress bar
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
