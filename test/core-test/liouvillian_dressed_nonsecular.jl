@testitem "Dressed Nonsecular Master Equation" begin
    using LinearAlgebra
    using SparseArrays

    function dressed_liouvillian(H, fields; tol = 1.0e-12)
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
    tol = 1.0e-6

    a = kron(destroy(N_c), qeye(2))
    sm = kron(qeye(N_c), sigmam())
    sp = sm'
    sx = sm + sp
    sz = sp * sm - sm * sp

    H = 1 * a' * a + 1 * sz / 2 + 0.5 * (a + a') * sx

    fields = [sqrt(0.01) * (a + a'), sqrt(0.01) * sx]
    Tlist = [0, 0.0]

    E, U, L1 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol)
    Ω = to_sparse((E' .- E)[1:N_trunc, 1:N_trunc], tol)

    H_d = U' * H * U
    Xp = QuantumObject(Ω .* triu(U' * (a + a') * U, 1).data, type = Operator())
    a2 = U' * a * U
    sm2 = U' * sm * U

    # Standard liouvillian case
    c_ops = [sqrt(0.01) * a2, sqrt(0.01) * sm2]
    L2 = liouvillian(H_d, c_ops)

    @test (expect(Hermitian(Xp' * Xp), steadystate(L1)) < 1.0e-10 && expect(Hermitian(Xp' * Xp), steadystate(L2)) > 1.0e-3)

    # Test that in the limit of σ_filter → 0, we recover the dressed master equation as in https://doi.org/10.1103/PhysRevA.84.043832
    L_gme_0 = liouvillian_dressed_nonsecular(H, fields, Tlist, tol = tol, σ_filter = 1.0e-14)[3]
    L_dr = dressed_liouvillian(H, fields, tol = tol)
    @test norm(L_gme_0.data - L_dr.data) / size(L_gme_0.data, 1) < 1.0e-8

    H = 1 * a' * a + 1 * sz / 2 + 1.0e-5 * (a * sp + a' * sm)

    Tlist = [0.2, 0.0]

    E, U, L1 = liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol)
    Ω = to_sparse((E' .- E)[1:N_trunc, 1:N_trunc], tol)

    H_d = U' * H * U
    Xp = QuantumObject(Ω .* triu(U' * (a + a') * U, 1).data, type = Operator())
    a2 = U' * a * U
    sm2 = U' * sm * U

    @test abs(expect(Xp' * Xp, steadystate(L1)) - n_thermal(1, Tlist[1])) / n_thermal(1, Tlist[1]) < 1.0e-4

    if VERSION >= v"1.11" # eigen is type unstable on v1.10
        @testset "Type Inference (liouvillian_dressed_nonsecular)" begin
            N_c = 30
            N_trunc = 10
            tol = 1.0e-14

            a = kron(destroy(N_c), qeye(2))
            sm = kron(qeye(N_c), sigmam())
            sp = sm'
            sx = sm + sp
            sz = sp * sm - sm * sp

            H = 1 * a' * a + 1 * sz / 2 + 0.5 * (a + a') * sx

            fields = [sqrt(0.01) * (a + a'), sqrt(0.01) * sx]
            Tlist = [0, 0.01]

            @inferred liouvillian_dressed_nonsecular(H, fields, Tlist, N_trunc = N_trunc, tol = tol)
        end
    end
end
