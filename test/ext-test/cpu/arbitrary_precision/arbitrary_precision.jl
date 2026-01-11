@testset "Arbitrary Precision" verbose = true begin
    N = 20
    Δ = 1.0
    U = 0.1
    F = 0.1
    κ = 0.1

    a = destroy(N)
    ψ0 = fock(N, 1)
    H = Δ * a' * a + U / 2 * a' * a' * a * a + F * (a + a')

    c_ops = [sqrt(κ) * a]

    tlist = range(0, 40, 100)

    a_big = SparseMatrixCSC{Complex{BigFloat}}(a)
    ψ0_big = Vector{Complex{BigFloat}}(ψ0)
    H_big = Δ * a_big' * a_big + U / 2 * a_big' * a_big' * a_big * a_big + F * (a_big + a_big')

    c_ops_big = [sqrt(BigFloat(κ)) * a_big]

    @testset "Arbitrary Precision (sesolve)" begin
        sol = sesolve(H, ψ0, tlist; progress_bar = Val(false))
        sol_big = sesolve(H_big, ψ0_big, tlist; progress_bar = Val(false))

        @test eltype(sol_big.states[1]) == Complex{BigFloat}

        # Test all fidelities are close to 1
        @test all((x) -> isapprox(fidelity(x[1], x[2]), 1; atol = 1.0e-7), zip(sol.states, sol_big.states))

        @inferred sesolve(H_big, ψ0_big, tlist; progress_bar = Val(false))
    end

    @testset "Arbitrary Precision (mesolve)" begin
        sol = mesolve(H, ψ0, tlist, c_ops; progress_bar = Val(false))
        sol_big = mesolve(H_big, ψ0_big, tlist, c_ops_big; progress_bar = Val(false))

        @test eltype(sol_big.states[1]) == Complex{BigFloat}

        # Test all fidelities are close to 1
        @test all((x) -> isapprox(hilbert_dist(x[1], x[2]), 0; atol = 1.0e-10), zip(sol.states, sol_big.states))

        @inferred mesolve(H_big, ψ0_big, tlist, c_ops_big; progress_bar = Val(false))
    end

    @testset "Arbitrary Precision (eigenvalues)" begin
        L = liouvillian(H, c_ops)
        L_big = liouvillian(H_big, c_ops_big)

        vals, vecs = eigenstates(L; sparse = Val(true), sigma = 0.01, eigvals = 7, krylovdim = 30)
        vals_big, vecs_big = eigenstates(L_big; sparse = Val(true), sigma = 0.01, eigvals = 7, krylovdim = 30)

        vals_al, vecs_al = eigsolve_al(L, 1 \ (50 * κ), eigvals = 7, krylovdim = 30)
        vals_big_al, vecs_big_al = eigsolve_al(
            L_big,
            1 \ (50 * BigFloat(κ)),
            eigvals = 7,
            krylovdim = 30,
            maxiter = 3,
        )

        # Align eigenvalues
        idxs = [findmin(abs.(vals_big .- val))[2] for val in vals]
        idxs_al = [findmin(abs.(vals_big_al .- val))[2] for val in vals_al]

        @test eltype(vals_big) == Complex{BigFloat}

        @test vals ≈ vals_big[idxs] atol = 1.0e-7
        @test all(zip(vecs, vecs_big[idxs])) do (v, v_big)
            return isapprox(abs(dot(v.data, v_big.data)), 1; atol = 1.0e-7)
        end

        @test vals_al[1:(end - 1)] ≈ vals_big_al[idxs_al][1:(end - 1)] atol = 1.0e-6
        @test all(zip(vecs_al, vecs_big_al[idxs_al])) do (v, v_big)
            return isapprox(abs(dot(v.data, v_big.data)), 1; atol = 1.0e-6)
        end
    end
end
