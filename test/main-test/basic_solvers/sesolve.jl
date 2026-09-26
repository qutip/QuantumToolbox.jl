using Test
using QuantumToolbox
import SciMLOperators: ScaledOperator

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "sesolve" begin

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

    @test prob.prob.f.f isa ScaledOperator
    @test !haskey(prob.prob.kwargs, :tstops) # tstops should not exist for time-independent cases
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
    @test sol.expect[1, saveat_idxs] ≈ expect(e_ops[1], sol3.states) atol = 1.0e-6

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
        @test allocs_tot < 95
    end

    @testset "Type Inference sesolve" begin
        @inferred sesolveProblem(H, ψ0, tlist, progress_bar = Val(false))
        @inferred sesolveProblem(H, ψ0, [0, 10], progress_bar = Val(false))
        @inferred sesolveProblem(H, TESetup.ψ0_int, tlist, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, progress_bar = Val(true)) # test progress bar
        @inferred sesolve(H, ψ0, tlist, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = e_ops, saveat = saveat, progress_bar = Val(false))
        @inferred sesolve(H, ψ0, tlist, e_ops = (TESetup.a' * TESetup.a, TESetup.a'), progress_bar = Val(false)) # We test the type inference for Tuple of different types
    end
end

@testset "sesolve_map" begin

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

    sols0 = sesolve_map(TESetup.H, ψ0_list, tlist; e_ops = e_ops, progress_bar = Val(false)) # no params
    sols1 = sesolve_map(H, ψ_0_e, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    sols2 = sesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    @test size(sols0) == (2,)
    @test sols0 isa Vector{<:TimeEvolutionSol}
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

            @test sol_0_e.expect[1, :] ≈ amp_rabi .* sin.(Ω_rabi * tlist) .^ 2 atol = 1.0e-2
            @test sol_1_g.expect[1, :] ≈ 1 .- amp_rabi .* sin.(Ω_rabi * tlist) .^ 2 atol = 1.0e-2
        end
    end

    @testset "Type Inference sesolve_map" begin
        @inferred sesolve_map(TESetup.H, ψ0_list, tlist; e_ops = e_ops, progress_bar = Val(true)) # no params, test progress bar
        @inferred sesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    end
end
