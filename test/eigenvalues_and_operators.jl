@testset "Eigenvalues and Operators" begin
    σx = sigmax()
    result = eigenstates(σx, sparse = false)
    λd, ψd, Td = result
    resstring = sprint((t, s) -> show(t, "text/plain", s), result)
    valstring = sprint((t, s) -> show(t, "text/plain", s), result.values)
    vecsstring = sprint((t, s) -> show(t, "text/plain", s), result.vectors)
    λs, ψs, Ts = eigenstates(σx, sparse = true, k = 2)
    λs1, ψs1, Ts1 = eigenstates(σx, sparse = true, k = 1)

    @test all([ψ.type isa KetQuantumObject for ψ in ψd])
    @test typeof(Td) <: AbstractMatrix
    @test typeof(Ts) <: AbstractMatrix
    @test typeof(Ts1) <: AbstractMatrix
    @test all(abs.(eigenenergies(σx, sparse = false)) .≈ abs.(λd))
    @test all(abs.(eigenenergies(σx, sparse = true, k = 2)) .≈ abs.(λs))
    @test resstring ==
          "EigsolveResult:   type=$(Operator)   dims=$(result.dims)\nvalues:\n$(valstring)\nvectors:\n$vecsstring"

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

    vals_d, vecs_d, mat_d = eigenstates(H_d)
    vals_c, vecs_c, mat_c = eigenstates(H_c)
    vals2, vecs2, mat2 = eigenstates(H_d, sparse = true, sigma = -0.9, k = 10, krylovdim = 30)
    sort!(vals_c, by = real)
    sort!(vals2, by = real)

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
    c_ops = [√((1 + n_thermal) * κ) * a, √κ * b, √(n_thermal * κ) * a_d]
    L = liouvillian(H, c_ops)

    # eigen solve for general matrices
    vals, _, vecs = eigsolve(L.data, sigma = 0.01, k = 10, krylovdim = 50)
    vals2, vecs2 = eigen(sparse_to_dense(L.data))
    vals3, state3, vecs3 = eigsolve_al(liouvillian(H, c_ops), 1 \ (40 * κ), k = 10, krylovdim = 50)
    idxs = sortperm(vals2, by = abs)
    vals2 = vals2[idxs][1:10]
    vecs2 = vecs2[:, idxs][:, 1:10]

    @test isapprox(sum(abs2, vals), sum(abs2, vals2), atol = 1e-7)
    @test isapprox(abs2(vals2[1]), abs2(vals3[1]), atol = 1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im * angle(vecs[1, 1])), vec2mat(vecs2[:, 1]), atol = 1e-7)
    @test isapprox(vec2mat(vecs[:, 1]) * exp(-1im * angle(vecs[1, 1])), vec2mat(vecs3[:, 1]), atol = 1e-5)

    # eigen solve for QuantumObject
    result = eigenstates(L, sparse = true, sigma = 0.01, k = 10, krylovdim = 50)
    vals, vecs = result
    resstring = sprint((t, s) -> show(t, "text/plain", s), result)
    valstring = sprint((t, s) -> show(t, "text/plain", s), result.values)
    vecsstring = sprint((t, s) -> show(t, "text/plain", s), result.vectors)
    @test resstring ==
          "EigsolveResult:   type=$(SuperOperator)   dims=$(result.dims)\nvalues:\n$(valstring)\nvectors:\n$vecsstring"

    vals2, vecs2 = eigenstates(L, sparse = false)
    idxs = sortperm(vals2, by = abs)
    vals2 = vals2[idxs][1:10]
    vecs2 = vecs2[idxs][1:10]

    @test result.type isa SuperOperatorQuantumObject
    @test result.dims == L.dims
    @test all([v.type isa OperatorKetQuantumObject for v in vecs])
    @test typeof(result.vectors) <: AbstractMatrix
    @test isapprox(sum(abs2, vals), sum(abs2, vals2), atol = 1e-7)
    @test isapprox(abs2(vals2[1]), abs2(vals3[1]), atol = 1e-7)
    @test isapprox(vec2mat(vecs[1]).data * exp(-1im * angle(vecs[1][1])), vec2mat(vecs2[1]).data, atol = 1e-7)
    @test isapprox(vec2mat(vecs[1]).data * exp(-1im * angle(vecs[1][1])), vec2mat(state3[1]).data, atol = 1e-5)
end
