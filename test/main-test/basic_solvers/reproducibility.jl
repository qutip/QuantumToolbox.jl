using Test
using QuantumToolbox
import Random: MersenneTwister

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "Reproducibility (mcsolve/ssesolve/smesolve)" begin

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

    @test sol_mc1.expect ≈ sol_mc2.expect atol = 1.0e-10
    @test sol_mc1.col_times ≈ sol_mc2.col_times atol = 1.0e-10
    @test sol_mc1.col_which ≈ sol_mc2.col_which atol = 1.0e-10

    @test sol_mc1.expect ≈ sol_mc3.expect[:, 1:500, :] atol = 1.0e-10

    @test sol_sse1.expect ≈ sol_sse2.expect atol = 1.0e-10

    @test sol_sse1.expect ≈ sol_sse3.expect[:, 1:50, :] atol = 1.0e-10

    @test sol_sme1.expect ≈ sol_sme2.expect atol = 1.0e-10

    @test sol_sme1.expect ≈ sol_sme3.expect[:, 1:50, :] atol = 1.0e-10
end