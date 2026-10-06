#!PTR_MULTITHREAD
using Test
using QuantumToolbox
import SciMLBase: EnsembleDistributed, EnsembleSplitThreads, EnsembleSerial

@testset "Distributed ensemble algorithms" begin
    @test Base.Threads.nthreads() > 1 # make sure the test runs with multi-threading

    # even though we are not executing the tests in parallel (distributed computing),
    # this test still checks the following:
    #   - correct display of progress bar across remote channels
    #   - correct output of the distributed solvers
    #
    # Note that we cannot use `@inferred` here: `solve` from SciMLBase relies on `pmap` for these ensemble algorithms,
    # and its return type is not inferrable (`Any`), even though the actual output is concrete.
    N = 10
    a = destroy(N)

    H = a' * a + 0.1 * (a^2 + a'^2)
    c_ops = (sqrt(0.05) * a,)
    sc_ops = c_ops

    L = liouvillian(H, c_ops)

    ψ0 = rand_ket(N)
    e_ops = (a' * a,)

    tlist = range(0, 100, 100)
    ntraj = 5

    # the deterministic solvers must return the same results independently of the ensemble algorithm
    sol_se_ref = sesolve_map(H, ψ0, tlist; e_ops, ensemblealg = EnsembleSerial(), progress_bar = Val(false))
    sol_me_ref = mesolve_map(L, ψ0, tlist; e_ops, ensemblealg = EnsembleSerial(), progress_bar = Val(false))

    for progress_bar in (Val(false), Val(true))
        for ensemblealg in (EnsembleDistributed(), EnsembleSplitThreads())
            sol_mc = mcsolve(H, ψ0, tlist, c_ops; ntraj, e_ops, ensemblealg, progress_bar)
            sol_sse = ssesolve(H, ψ0, tlist, c_ops; ntraj, e_ops, ensemblealg, progress_bar)
            sol_sme = smesolve(H, ψ0, tlist, c_ops, sc_ops; ntraj, e_ops, ensemblealg, progress_bar)

            @test sol_mc isa TimeEvolutionMCSol
            @test sol_sse isa TimeEvolutionStochasticSol
            @test sol_sme isa TimeEvolutionStochasticSol
            for sol in (sol_mc, sol_sse, sol_sme)
                @test sol.ntraj == ntraj
                @test size(sol.expect) == (length(e_ops), length(tlist))
                @test all(isfinite, sol.expect)
            end

            sol_se = sesolve_map(H, ψ0, tlist; e_ops, ensemblealg, progress_bar)
            sol_me = mesolve_map(L, ψ0, tlist; e_ops, ensemblealg, progress_bar)
            sol_mc_map = mcsolve_map(H, ψ0, tlist, c_ops; ntraj, e_ops, ensemblealg, progress_bar)

            @test size(sol_se) == size(sol_me) == size(sol_mc_map) == (1,)
            @test sol_se[1].expect ≈ sol_se_ref[1].expect
            @test sol_me[1].expect ≈ sol_me_ref[1].expect
            @test sol_mc_map[1] isa TimeEvolutionMCSol
            @test sol_mc_map[1].ntraj == ntraj
            @test size(sol_mc_map[1].expect) == (length(e_ops), length(tlist))
            @test all(isfinite, sol_mc_map[1].expect)
        end
    end
end
