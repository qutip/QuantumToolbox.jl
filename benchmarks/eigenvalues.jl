function benchmark_eigenvalues()
    N = 5
    a = tensor(destroy(N), qeye(N))
    a_d = a'
    b = tensor(qeye(N), destroy(N))
    b_d = b'

    ωc = 1
    ωb = 1
    g = 0.2
    κ = 0.01
    n_thermal = 0.1

    H = ωc * a_d * a + ωb * b_d * b + g * (a + a_d) * (b + b_d)
    c_ops = [√((1 + n_thermal) * κ) * a, √κ * b, √(n_thermal * κ) * a_d]
    L = liouvillian(H, c_ops)

    SUITE["Eigenvalues"]["eigenstates"]["dense"] = @benchmarkable eigenstates($L)
    return SUITE["Eigenvalues"]["eigenstates"]["sparse"] =
        @benchmarkable eigenstates($L, sparse = true, sigma = 0.01, k = 5)
end

benchmark_eigenvalues()
