function benchmark_timeevolution()
    ωc = 1
    ωq = 1
    g = 0.1
    ωd = 0.99
    F = 0.07

    Δc = ωc - ωd
    Δq = ωq - ωd

    # Operators definition
    N = 20 # cutoff for the cavity Hilbert space
    a = tensor(destroy(N), qeye(2))
    σm = tensor(qeye(N), sigmam())
    σz = tensor(qeye(N), sigmaz())

    # Hamiltonian
    H = Δc * a' * a + Δq * σz / 2 + g * (a' * σm + a * σm') + F * (a + a')

    e_ops = [a' * a, σz]

    # Initial state
    ψ0 = tensor(coherent(N, 0), fock(2, 1))

    ## sesolve ##

    tlist = range(0, 2π * 10 / g, 1000)

    SUITE["Time Evolution"]["time-independent"]["sesolve"] =
        @benchmarkable sesolve($H, $ψ0, $tlist, e_ops = $e_ops, progress_bar = false)

    ## mesolve ##

    nth = 0.01
    γ = 0.05
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a', sqrt(γ) * σm]

    tlist = range(0, 10 / γ, 100)

    SUITE["Time Evolution"]["time-independent"]["mesolve"] =
        @benchmarkable mesolve($H, $ψ0, $tlist, $c_ops, e_ops = $e_ops, progress_bar = false)

    ## mcsolve ##

    SUITE["Time Evolution"]["time-independent"]["mcsolve"]["Serial"] = @benchmarkable mcsolve(
        $H,
        $ψ0,
        $tlist,
        $c_ops,
        n_traj = 100,
        e_ops = $e_ops,
        progress_bar = false,
        ensemble_method = EnsembleSerial(),
    )
    return SUITE["Time Evolution"]["time-independent"]["mcsolve"]["Multithreaded"] = @benchmarkable mcsolve(
        $H,
        $ψ0,
        $tlist,
        $c_ops,
        n_traj = 100,
        e_ops = $e_ops,
        progress_bar = false,
        ensemble_method = EnsembleThreads(),
    )
end

benchmark_timeevolution()
