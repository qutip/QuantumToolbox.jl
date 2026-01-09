function benchmark_dsf!(SUITE)
    F = 2
    Δ = 0.25
    κ = 1
    U = 0.01
    J = 0.5
    tlist = LinRange(0, 15, 300)

    N = 5
    a1 = kron(destroy(N), qeye(N))
    a2 = kron(qeye(N), destroy(N))
    function H_dsf(op_list, p)
        Δ = p.Δ
        F = p.F
        U = p.U
        J = p.J
        a1, a2 = op_list
        return Δ * a1' * a1 +
            Δ * a2' * a2 +
            U * a1'^2 * a1^2 +
            U * a2'^2 * a2^2 +
            F * (a1 + a1') +
            J * (a1' * a2 + a1 * a2')
    end
    function c_ops_dsf(op_list, p)
        κ = p.κ
        a1, a2 = op_list
        return [√κ * a1, √κ * a2]
    end
    function e_ops_dsf(op_list, p)
        a1, a2 = op_list
        return [a1' * a1, a2' * a2]
    end
    op_list = [a1, a2]
    ψ0 = kron(fock(N, 0), fock(N, 0))
    α0 = 1.5
    α0_l = [α0, α0]
    dsf_params = (Δ = Δ, F = F, κ = κ, U = U, J = J)

    SUITE["Time Evolution"]["Dynamical Shifted Fock"]["mesolve"] = @benchmarkable dsf_mesolve(
        $H_dsf,
        $ψ0,
        $tlist,
        $c_ops_dsf,
        $op_list,
        $α0_l,
        $dsf_params,
        e_ops = $e_ops_dsf,
        progress_bar = Val(false),
    )

    SUITE["Time Evolution"]["Dynamical Shifted Fock"]["mcsolve"]["Serial"] = @benchmarkable dsf_mcsolve(
        $H_dsf,
        $ψ0,
        $tlist,
        $c_ops_dsf,
        $op_list,
        $α0_l,
        $dsf_params,
        ntraj = 100,
        e_ops = $e_ops_dsf,
        progress_bar = Val(false),
        ensemblealg = EnsembleSerial(),
    )

    SUITE["Time Evolution"]["Dynamical Shifted Fock"]["mcsolve"]["Multithreaded"] = @benchmarkable dsf_mcsolve(
        $H_dsf,
        $ψ0,
        $tlist,
        $c_ops_dsf,
        $op_list,
        $α0_l,
        $dsf_params,
        ntraj = 100,
        e_ops = $e_ops_dsf,
        progress_bar = Val(false),
        ensemblealg = EnsembleThreads(),
    )

    return nothing
end
