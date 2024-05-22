function H_dfd2(dims, p)
    Δ = p.Δ
    F = p.F
    J = p.J
    a = tensor(destroy(dims[1]), qeye(dims[2]))
    b = tensor(qeye(dims[1]), destroy(dims[2]))
    return Δ * a' * a + F * (a + a') + Δ * b' * b + J * (a' * b + a * b')
end
function c_ops_dfd2(dims, p)
    κ = p.κ
    a = tensor(destroy(dims[1]), qeye(dims[2]))
    b = tensor(qeye(dims[1]), destroy(dims[2]))
    return [√κ * a, √κ * b]
end
function e_ops_dfd2(dims, p)
    a = tensor(destroy(dims[1]), qeye(dims[2]))
    b = tensor(qeye(dims[1]), destroy(dims[2]))
    return [a' * a, b' * b]
end

function benchmark_dfd()
    F, Δ, κ, J = 1, 0.25, 1, 0.05
    maxdims = [50, 50]

    ψ0 = tensor(fock(3, 0), fock(20, 15))
    dfd_params = (Δ = Δ, F = F, κ = κ, J = J)

    tlist = range(0, 15 / κ, 100)

    return SUITE["Time Evolution"]["Dynamical Fock Dimension"] = @benchmarkable dfd_mesolve(
        H_dfd2,
        $ψ0,
        $tlist,
        c_ops_dfd2,
        $maxdims,
        $dfd_params,
        e_ops = e_ops_dfd2,
        progress_bar = false,
    )
end

benchmark_dfd()
