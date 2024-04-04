@testset "Eigenvalues and Operators" begin
    N = 30
    a = kron(destroy(N), qeye(2))
    a_d = a'

    sm = kron(qeye(N), sigmam())
    sp = sm'
    sx = kron(qeye(N), sigmax())
    sy = kron(qeye(N), sigmay())
    sz = kron(qeye(N), sigmaz())

    η = 0.2
    H_d = a_d * a + 0.5 * sz - 1im * η * (a - a_d) * sx + η^2
    H_c = a_d * a + 0.5 * (sz * cosm(2 * η * (a + a_d)) + sy * sinm(2 * η * (a + a_d)))

    vals_d, vecs_d = eigen(H_d)
    vals_c, vecs_c = eigen((H_c))
    vals2, vecs2 = eigsolve(H_d, sigma=-0.9, k=10, krylovdim=30)
    sort!(vals_c, by=real)
    sort!(vals2, by=real)

    @test sum(real.(vals_d[1:20]) .- real.(vals_c[1:20])) / 20 < 1e-3
    @test sum(real.(vals_d[1:10]) .- real.(vals2[1:10])) / 20 < 1e-3


    N = 5
    a = kron(destroy(N), qeye(N))
    a_d = a'
    b = kron(qeye(N), destroy(N))
    b_d = b'

    ωc = 1
    ωb = 1
    g = 0.01
    κ = 0.1
    n_thermal = 0.01

    H = ωc * a_d * a + ωb * b_d * b + g * (a + a_d) * (b + b_d)
    c_ops = [√((1+n_thermal)*κ) * a, √κ * b, √(n_thermal*κ) * a_d]
    L = liouvillian(H, c_ops).data

    vals, vecs = eigsolve(L, sigma=0.01, k=10, krylovdim=50)
    vals2, vecs2 = eigen(sparse_to_dense(L))
    vals3, vecs3 = eigsolve_al(liouvillian(H, c_ops), 1\(40*κ), k=10, krylovdim=50)
    idxs = sortperm(vals2, by=abs)
    vals2 = vals2[idxs][1:10]
    vecs2 = vecs2[:, idxs][:, 1:10]

    @test isapprox(sum(abs2, vals), sum(abs2, vals2), atol=1e-7)
    @test isapprox(abs2(vals2[1]), abs2(vals3[1]), atol=1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im*angle(vecs[1,1])), vec2mat(vecs2[:, 1]), atol=1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im*angle(vecs[1,1])), vec2mat(vecs3[:, 1]), atol=1e-5)
end