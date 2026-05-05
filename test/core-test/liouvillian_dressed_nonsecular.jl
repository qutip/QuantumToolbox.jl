@testitem "Dressed Nonsecular Master Equation" begin
    using LinearAlgebra
    using SparseArrays

    function dressed_liouvillian(H, fields; tol = 1e-12)
        vals, vecs = eigenstates(H)
        N_tot = size(H, 1)

        H_d = QuantumObject(spdiagm(0 => complex.(vals)), type = Operator())

        L = liouvillian(H_d)

        iter = Iterators.product(1:N_tot, 1:N_tot)
        iter_filtered = Iterators.filter(x -> x[2] > x[1], iter)
        foreach(iter_filtered) do idx
            i, j = idx
            Γ_ij = sum(eachindex(fields)) do k
                (vals[j] - vals[i]) * abs2(matrix_element(vecs[i], fields[k], vecs[j]))
            end
            if Γ_ij > tol
                L += Γ_ij * lindblad_dissipator(projection(N_tot, i - 1, j - 1))
            end
        end

        return L
    end

    N_c = 30
    N_trunc = 10
    tol = 1e-6

    a = kron(destroy(N_c), qeye(2))
    sm = kron(qeye(N_c), sigmam())
    sp = sm'
    sx = sm + sp
    sz = sp * sm - sm * sp

    H = 1 * a' * a + 1 * sz / 2 + 0.5 * im * (a - a') * sx

    fields = (sqrt(0.01) * im * (a - a'), sqrt(0.01) * sx)
    Tlist = (0.0, 0.0)
    σ_filter = 4.0

    E, U, L1 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = σ_filter)
    Ω = to_sparse((E' .- E)[1:N_trunc, 1:N_trunc], tol)

    H_d = U' * H * U
    Xp = QuantumObject(Ω .* triu(U' * (a + a') * U, 1).data, type = Operator())
    a2 = U' * a * U
    sm2 = U' * sm * U

    # Standard liouvillian case
    c_ops = [sqrt(0.01) * a2, sqrt(0.01) * sm2]
    L_std = liouvillian(H_d, c_ops)

    @test (expect(Hermitian(Xp' * Xp), steadystate(L1)) < 1e-10 && expect(Hermitian(Xp' * Xp), steadystate(L_std)) > 1e-3)

    # Test that in the limit of σ_filter → 0, we recover the dressed master equation as in https://doi.org/10.1103/PhysRevA.84.043832
    L_gme_0 = liouvillian_dressed_nonsecular(H, fields, Tlist, tol = tol, σ_filter = 1e-14)[3]
    L_dr = dressed_liouvillian(H, fields, tol = tol)
    @test norm(L_gme_0.data - L_dr.data) / size(L_gme_0.data, 1) < 1.0e-8

    # Test that the σ = ∞ and σ >> 1 limits give the same result
    L1 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = 1.0e5)[3]
    L2 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = Inf)[3]
    @test norm(L1.data - L2.data) / size(L1.data, 1) < 1.0e-7

    # Test matrix-form output against vectorized output for the unfiltered case.
    L_vec_nf = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = nothing)[3]
    L_mat_nf =
        liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = nothing, matrix_form = Val(true))[3]

    ρ_test = rand_dm(N_trunc)
    dρ_vec = vector_to_operator(L_vec_nf * operator_to_vector(ρ_test))
    L_mat_nf_cached = cache_operator(L_mat_nf, ρ_test)
    dρ_mat = L_mat_nf_cached * ρ_test
    @test dρ_mat ≈ dρ_vec atol = 1e-10

    H = 1 * a' * a + 1 * sz / 2 + 1e-5 * (a * sp + a' * sm)

    Tlist = (0.2, 0.0)

    E, U, L1 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = σ_filter)
    Ω = to_sparse((E' .- E)[1:N_trunc, 1:N_trunc], tol)

    H_d = U' * H * U
    Xp = QuantumObject(Ω .* triu(U' * (a + a') * U, 1).data, type = Operator())
    a2 = U' * a * U
    sm2 = U' * sm * U

    @test abs(expect(Xp' * Xp, steadystate(L1)) - n_thermal(1, Tlist[1])) / n_thermal(1, Tlist[1]) < 1e-4

    if VERSION >= v"1.11" # eigen is type unstable on v1.10
        @testset "Type Inference (liouvillian_dressed_nonsecular)" begin
            N_c = 30
            N_trunc = 10
            tol = 1e-14

            a = kron(destroy(N_c), qeye(2))
            sm = kron(qeye(N_c), sigmam())
            sp = sm'
            sx = sm + sp
            sz = sp * sm - sm * sp

            H = 1 * a' * a + 1 * sz / 2 + 0.5 * (a + a') * sx

            fields = (sqrt(0.01) * im * (a - a'), sqrt(0.01) * sx)
            Tlist = (0.0, 0.01)

            σ_filter = 4.0

            @inferred liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = nothing)
            @inferred liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = σ_filter)

            @inferred liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol, σ_filter = nothing, matrix_form = Val(true))
        end
    end
end
