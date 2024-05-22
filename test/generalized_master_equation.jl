@testset "Generalized Master Equation" begin
    N_c = 30
    N_trunc = 10
    tol = 1e-14

    a = kron(destroy(N_c), qeye(2))
    sm = kron(qeye(N_c), sigmam())
    sp = sm'
    sx = sm + sp
    sz = sp * sm - sm * sp

    H = 1 * a' * a + 1 * sz / 2 + 0.5 * (a + a') * sx

    fields = [sqrt(0.01) * (a + a'), sqrt(0.01) * sx]
    Tlist = [0, 0.0]

    E, U, L1 = liouvillian_generalized(H, fields, Tlist, N_trunc = N_trunc, tol = tol)
    Ω = dense_to_sparse((E'.-E)[1:N_trunc, 1:N_trunc], tol)

    H_d = Qobj(dense_to_sparse((U'*H*U)[1:N_trunc, 1:N_trunc], tol))
    Xp = Qobj(Ω .* dense_to_sparse(triu((U'*(a+a')*U).data[1:N_trunc, 1:N_trunc], 1), tol))
    a2 = Qobj(dense_to_sparse((U'*a*U).data[1:N_trunc, 1:N_trunc], tol))
    sm2 = Qobj(dense_to_sparse((U'*sm*U).data[1:N_trunc, 1:N_trunc], tol))

    # Standard liouvillian case
    c_ops = [sqrt(0.01) * a2, sqrt(0.01) * sm2]
    L2 = liouvillian(H_d, c_ops)

    @test (expect(Xp' * Xp, steadystate(L1)) < 1e-10 && expect(Xp' * Xp, steadystate(L2)) > 1e-3)

    H = 1 * a' * a + 1 * sz / 2 + 1e-5 * (a * sp + a' * sm)

    Tlist = [0.2, 0.0]

    E, U, L1 = liouvillian_generalized(H, fields, Tlist, N_trunc = N_trunc, tol = tol)
    Ω = dense_to_sparse((E'.-E)[1:N_trunc, 1:N_trunc], tol)

    H_d = Qobj(dense_to_sparse((U'*H*U)[1:N_trunc, 1:N_trunc], tol))
    Xp = Qobj(Ω .* dense_to_sparse(triu((U'*(a+a')*U).data[1:N_trunc, 1:N_trunc], 1), tol))
    a2 = Qobj(dense_to_sparse((U'*a*U).data[1:N_trunc, 1:N_trunc], tol))
    sm2 = Qobj(dense_to_sparse((U'*sm*U).data[1:N_trunc, 1:N_trunc], tol))

    @test abs(expect(Xp' * Xp, steadystate(L1)) - n_th(1, Tlist[1])) / n_th(1, Tlist[1]) < 1e-4
end
