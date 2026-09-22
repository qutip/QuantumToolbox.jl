using Test
using QuantumToolbox
import SciMLOperators: MatrixOperator

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "mesolve" begin

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
    prob_me_mat = mesolveProblem(H, ket2dm(ψ0), tlist, c_ops, progress_bar = Val(false), matrix_form = Val(true))
    sol_me_mat = mesolve(H, ket2dm(ψ0), tlist, c_ops, progress_bar = Val(false), matrix_form = Val(true))
    sol_me_mat2 =
        mesolve(H, ket2dm(ψ0), tlist, c_ops, e_ops = e_ops, saveat = saveat, progress_bar = Val(false), matrix_form = Val(true))

    # For testing the `OperatorKet` input
    sol_me4 = mesolve(H, operator_to_vector(ket2dm(ψ0)), tlist, c_ops, saveat = saveat, progress_bar = Val(false))

    # Redirect to `sesolve`
    sol_me5 = mesolve(H, ψ0, tlist, progress_bar = Val(false))

    @test TESetup.prob_me.prob.f.f isa MatrixOperator
    @test !haskey(TESetup.prob_me.prob.kwargs, :tstops) # tstops should not exist for time-independent cases
    @test !haskey(prob_me_mat.prob.kwargs, :tstops)
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
    @test sol_me3.expect[1, TESetup.saveat_idxs] ≈ expect(e_ops[1], sol_me3.states) atol = 1.0e-6
    @test all([sol_me3.states[i] ≈ vector_to_operator(sol_me4.states[i]) for i in eachindex(saveat)])
    @test length(sol_me_mat.times_states) == length(tlist)
    @test length(sol_me_mat.states) == length(tlist)
    @test sol_me_mat.expect === nothing
    @test all(isoper, sol_me_mat.states)
    @test size(sol_me_mat2.expect) == (length(e_ops), length(tlist))
    @test sol_me_mat2.expect[1, TESetup.saveat_idxs] ≈ expect(e_ops[1], sol_me_mat2.states) atol = 1.0e-6
    @test sol_me_mat2.expect[1, :] ≈ sol_me.expect[1, :] atol = 1.0e-6

    @testset "Pure-dissipator Liouvillian matrix form" begin
        L_diss_mat = liouvillian(nothing, c_ops; matrix_form = Val(true))
        L_diss_mat_vector = liouvillian(nothing, collect(c_ops); matrix_form = Val(true))
        @test L_diss_mat.type isa SuperOperatorMatrixForm
        @test L_diss_mat_vector.type isa SuperOperatorMatrixForm
    end

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
        @inferred mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), matrix_form = Val(true))
        @inferred liouvillian(nothing, c_ops; matrix_form = Val(true))
        @inferred mesolveProblem(H, ψ0, [0, 10], c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolveProblem(H, TESetup.ψ0_int, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(true)) # also test progress bar
        @inferred mesolve(H, ψ0, tlist, c_ops, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), matrix_form = Val(true))
        @inferred mesolve(H, ψ0, tlist, c_ops, e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
        @inferred mesolve(H, ψ0, tlist, (a, ad_t), e_ops = (a' * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple
        @inferred mesolve(TESetup.H_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        @inferred mesolve(TESetup.H_td2, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
        @inferred mesolve(TESetup.L_td, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false), params = p)
    end
end

@testset "mesolve_map" begin

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

    # Test with multiple initial states but no params
    sols0 = mesolve_map(TESetup.H, ψ0_list, tlist, c_ops; e_ops = e_ops, progress_bar = Val(false))
    # Test with single initial state
    sols1 = mesolve_map(H, ψ_0_e, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    # Test with multiple initial states
    sols2 = mesolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))
    # Test matrix_form = Val(true) case
    sols2_mat = mesolve_map(H, ψ0_list, tlist, c_ops; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false), matrix_form = Val(true))

    # Test redirect to sesolve_map when c_ops is nothing
    sols3 = mesolve_map(H, ψ0_list, tlist; e_ops = e_ops, params = (ωc_list, ωq_list), progress_bar = Val(false))

    @test size(sols0) == (2,)
    @test sols0 isa Vector{<:TimeEvolutionSol}
    @test size(sols1) == (1, 3, 4)
    @test sols1 isa Array{<:TimeEvolutionSol}
    @test size(sols2) == (2, 3, 4)
    @test sols2 isa Array{<:TimeEvolutionSol}
    @test size(sols2_mat) == (2, 3, 4)
    @test sols2_mat isa Array{<:TimeEvolutionSol}
    @test size(sols3) == (2, 3, 4)
    @test sols3 isa Array{<:TimeEvolutionSol}

    # Verify that solutions make physical sense
    for (i, ωc) in enumerate(ωc_list)
        for (j, ωq) in enumerate(ωq_list)
            sol_0_e = sols2[1, i, j]
            sol_1_g = sols2[2, i, j]
            sol_0_e_mat = sols2_mat[1, i, j]
            sol_1_g_mat = sols2_mat[2, i, j]

            # Check that expectation values are bounded and physical (take real part for physical observables)
            @test all(x -> real(x) >= -1.0e-4, sol_0_e.expect[1, :]) # a'a should be non-negative (with small tolerance)
            @test all(x -> real(x) >= -1.0e-4, sol_1_g.expect[1, :])
            @test all(x -> real(x) >= -1.0e-4, sol_0_e_mat.expect[1, :])
            @test all(x -> real(x) >= -1.0e-4, sol_1_g_mat.expect[1, :])
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
        @inferred mesolve_map(TESetup.H, ψ0_list, tlist, c_ops; e_ops = e_ops, progress_bar = Val(true)) # no params, but test progress bar
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