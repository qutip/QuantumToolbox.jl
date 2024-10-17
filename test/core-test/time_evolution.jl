@testset "Time Evolution and Partial Trace" verbose = true begin
    @testset "sesolve" begin
        N = 10
        a_d = kron(create(N), qeye(2))
        a = a_d'
        sx = kron(qeye(N), sigmax())
        sy = tensor(qeye(N), sigmay())
        sz = qeye(N) ⊗ sigmaz()
        η = 0.01
        H = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx
        psi0 = kron(fock(N, 0), fock(2, 0))
        t_l = LinRange(0, 1000, 1000)
        e_ops = [a_d * a]
        sol = sesolve(H, psi0, t_l, e_ops = e_ops, progress_bar = Val(false))
        sol2 = sesolve(H, psi0, t_l, progress_bar = Val(false))
        sol3 = sesolve(H, psi0, t_l, e_ops = e_ops, saveat = t_l, progress_bar = Val(false))
        sol_string = sprint((t, s) -> show(t, "text/plain", s), sol)
        @test sum(abs.(sol.expect[1, :] .- sin.(η * t_l) .^ 2)) / length(t_l) < 0.1
        @test length(sol.times) == length(t_l)
        @test length(sol.states) == 1
        @test size(sol.expect) == (length(e_ops), length(t_l))
        @test length(sol2.times) == length(t_l)
        @test length(sol2.states) == length(t_l)
        @test size(sol2.expect) == (0, length(t_l))
        @test length(sol3.times) == length(t_l)
        @test length(sol3.states) == length(t_l)
        @test size(sol3.expect) == (length(e_ops), length(t_l))
        @test sol_string ==
              "Solution of time evolution\n" *
              "(return code: $(sol.retcode))\n" *
              "--------------------------\n" *
              "num_states = $(length(sol.states))\n" *
              "num_expect = $(size(sol.expect, 1))\n" *
              "ODE alg.: $(sol.alg)\n" *
              "abstol = $(sol.abstol)\n" *
              "reltol = $(sol.reltol)\n"

        @testset "Type Inference sesolve" begin
            @inferred sesolveProblem(H, psi0, t_l, progress_bar = Val(false))
            @inferred sesolveProblem(H, psi0, [0, 10], progress_bar = Val(false))
            @inferred sesolveProblem(H, Qobj(zeros(Int64, N * 2); dims = (N, 2)), t_l, progress_bar = Val(false))
            @inferred sesolve(H, psi0, t_l, e_ops = e_ops, progress_bar = Val(false))
            @inferred sesolve(H, psi0, t_l, progress_bar = Val(false))
            @inferred sesolve(H, psi0, t_l, e_ops = e_ops, saveat = t_l, progress_bar = Val(false))
            @inferred sesolve(H, psi0, t_l, e_ops = (a_d * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple of different types
        end
    end

    @testset "mesolve, mcsolve, and ssesolve" begin
        N = 10
        a = destroy(N)
        a_d = a'
        H = a_d * a
        c_ops = [sqrt(0.1) * a]
        e_ops = [a_d * a]
        psi0 = basis(N, 3)
        t_l = LinRange(0, 100, 1000)
        sol_me = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
        sol_me2 = mesolve(H, psi0, t_l, c_ops, progress_bar = Val(false))
        sol_me3 = mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, saveat = t_l, progress_bar = Val(false))
        sol_mc = mcsolve(H, psi0, t_l, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false))
        sol_mc_states = mcsolve(H, psi0, t_l, c_ops, ntraj = 500, saveat = t_l, progress_bar = Val(false))
        sol_sse = ssesolve(H, psi0, t_l, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false))

        ρt_mc = [ket2dm.(normalize.(states)) for states in sol_mc_states.states]
        expect_mc_states = mapreduce(states -> expect.(Ref(e_ops[1]), states), hcat, ρt_mc)
        expect_mc_states_mean = sum(expect_mc_states, dims = 2) / size(expect_mc_states, 2)

        sol_me_string = sprint((t, s) -> show(t, "text/plain", s), sol_me)
        sol_mc_string = sprint((t, s) -> show(t, "text/plain", s), sol_mc)
        sol_sse_string = sprint((t, s) -> show(t, "text/plain", s), sol_sse)
        @test sum(abs.(sol_mc.expect .- sol_me.expect)) / length(t_l) < 0.1
        @test sum(abs.(vec(expect_mc_states_mean) .- vec(sol_me.expect))) / length(t_l) < 0.1
        @test sum(abs.(sol_sse.expect .- sol_me.expect)) / length(t_l) < 0.1
        @test length(sol_me.times) == length(t_l)
        @test length(sol_me.states) == 1
        @test size(sol_me.expect) == (length(e_ops), length(t_l))
        @test length(sol_me2.times) == length(t_l)
        @test length(sol_me2.states) == length(t_l)
        @test size(sol_me2.expect) == (0, length(t_l))
        @test length(sol_me3.times) == length(t_l)
        @test length(sol_me3.states) == length(t_l)
        @test size(sol_me3.expect) == (length(e_ops), length(t_l))
        @test length(sol_mc.times) == length(t_l)
        @test size(sol_mc.expect) == (length(e_ops), length(t_l))
        @test length(sol_mc_states.times) == length(t_l)
        @test size(sol_mc_states.expect) == (0, length(t_l))
        @test length(sol_sse.times) == length(t_l)
        @test size(sol_sse.expect) == (length(e_ops), length(t_l))
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
        @test sol_sse_string ==
              "Solution of quantum trajectories\n" *
              "(converged: $(sol_sse.converged))\n" *
              "--------------------------------\n" *
              "num_trajectories = $(sol_sse.ntraj)\n" *
              "num_states = $(length(sol_sse.states[1]))\n" *
              "num_expect = $(size(sol_sse.expect, 1))\n" *
              "SDE alg.: $(sol_sse.alg)\n" *
              "abstol = $(sol_sse.abstol)\n" *
              "reltol = $(sol_sse.reltol)\n"

        @testset "Type Inference mesolve" begin
            @inferred mesolveProblem(H, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolveProblem(H, psi0, [0, 10], c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolveProblem(H, Qobj(zeros(Int64, N)), t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, progress_bar = Val(false))
            @inferred mesolve(H, psi0, t_l, c_ops, progress_bar = Val(false))
            @inferred mesolve(H, psi0, t_l, c_ops, e_ops = e_ops, saveat = t_l, progress_bar = Val(false))
            @inferred mesolve(H, psi0, t_l, (a, a'), e_ops = (a_d * a, a'), progress_bar = Val(false)) # We test the type inference for Tuple of different types
        end

        @testset "Type Inference mcsolve" begin
            @inferred mcsolveEnsembleProblem(H, psi0, t_l, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false))
            @inferred mcsolve(H, psi0, t_l, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false))
            @inferred mcsolve(H, psi0, t_l, c_ops, ntraj = 500, progress_bar = Val(true))
            @inferred mcsolve(H, psi0, [0, 10], c_ops, ntraj = 500, progress_bar = Val(false))
            @inferred mcsolve(H, Qobj(zeros(Int64, N)), t_l, c_ops, ntraj = 500, progress_bar = Val(false))
            @inferred mcsolve(H, psi0, t_l, (a, a'), e_ops = (a_d * a, a'), ntraj = 500, progress_bar = Val(false)) # We test the type inference for Tuple of different types
        end

        @testset "Type Inference ssesolve" begin
            @inferred ssesolveEnsembleProblem(
                H,
                psi0,
                t_l,
                c_ops,
                ntraj = 500,
                e_ops = e_ops,
                progress_bar = Val(false),
            )
            @inferred ssesolve(H, psi0, t_l, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false))
            @inferred ssesolve(H, psi0, t_l, c_ops, ntraj = 500, progress_bar = Val(true))
        end

        @testset "mcsolve and ssesolve reproducibility" begin
            N = 10
            a = tensor(destroy(N), qeye(2))
            σm = tensor(qeye(N), sigmam())
            σp = σm'
            σz = tensor(qeye(N), sigmaz())

            ω = 1.0
            g = 0.1
            γ = 0.01
            nth = 0.1

            H = ω * a' * a + ω * σz / 2 + g * (a' * σm + a * σp)
            c_ops = [sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a', sqrt(γ * (1 + nth)) * σm, sqrt(γ * nth) * σp]
            e_ops = [a' * a, σz]

            psi0 = tensor(basis(N, 0), basis(2, 0))
            tlist = range(0, 20 / γ, 1000)

            rng = MersenneTwister(1234)
            sleep(0.1) # If we don't sleep, we get an error (why?)
            sol_mc1 = mcsolve(H, psi0, tlist, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false), rng = rng)
            sol_sse1 = ssesolve(H, psi0, tlist, c_ops, ntraj = 50, e_ops = e_ops, progress_bar = Val(false), rng = rng)

            rng = MersenneTwister(1234)
            sleep(0.1)
            sol_mc2 = mcsolve(H, psi0, tlist, c_ops, ntraj = 500, e_ops = e_ops, progress_bar = Val(false), rng = rng)
            sol_sse2 = ssesolve(H, psi0, tlist, c_ops, ntraj = 50, e_ops = e_ops, progress_bar = Val(false), rng = rng)

            rng = MersenneTwister(1234)
            sleep(0.1)
            sol_mc3 = mcsolve(H, psi0, tlist, c_ops, ntraj = 510, e_ops = e_ops, progress_bar = Val(false), rng = rng)

            @test sol_mc1.expect ≈ sol_mc2.expect atol = 1e-10
            @test sol_mc1.expect_all ≈ sol_mc2.expect_all atol = 1e-10
            @test sol_mc1.jump_times ≈ sol_mc2.jump_times atol = 1e-10
            @test sol_mc1.jump_which ≈ sol_mc2.jump_which atol = 1e-10

            @test sol_mc1.expect_all ≈ sol_mc3.expect_all[1:500, :, :] atol = 1e-10

            @test sol_sse1.expect ≈ sol_sse2.expect atol = 1e-10
            @test sol_sse1.expect_all ≈ sol_sse2.expect_all atol = 1e-10
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
        sol_mc = mcsolve(H, psi0, t_l, c_ops, ntraj = 500, e_ops = [sp1 * sm1, sp2 * sm2], progress_bar = Val(false))
        @test sum(abs.(sol_mc.expect[1:2, :] .- sol_me.expect[1:2, :])) / length(t_l) < 0.1
        @test expect(sp1 * sm1, sol_me.states[end]) ≈ expect(sigmap() * sigmam(), ptrace(sol_me.states[end], 1))
    end
end
