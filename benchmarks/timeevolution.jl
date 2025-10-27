function benchmark_timeevolution!(SUITE)
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
        @benchmarkable sesolve($H, $ψ0, $tlist, e_ops = $e_ops, progress_bar = Val(false))

    ## mesolve ##

    nth = 0.01
    γ = 0.05
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a', sqrt(γ) * σm]

    tlist = range(0, 10 / γ, 100)

    SUITE["Time Evolution"]["time-independent"]["mesolve"] =
        @benchmarkable mesolve($H, $ψ0, $tlist, $c_ops, e_ops = $e_ops, progress_bar = Val(false))

    ## mcsolve ##

    SUITE["Time Evolution"]["time-independent"]["mcsolve"]["Serial"] = @benchmarkable mcsolve(
        $H,
        $ψ0,
        $tlist,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        ensemblealg = EnsembleSerial(),
    )
    SUITE["Time Evolution"]["time-independent"]["mcsolve"]["Multithreaded"] = @benchmarkable mcsolve(
        $H,
        $ψ0,
        $tlist,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        ensemblealg = EnsembleThreads(),
    )

    ## Time-dependent evolutions ##

    # Hamiltonian in the lab frame (without drive frame transformation)
    H_lab = ωc * a' * a + ωq / 2 * σz + g * (a' * σm + a * σm')

    # Define time-dependent drive terms
    coef1(p, t) = p.F * exp(1im * p.ωd * t)
    coef2(p, t) = p.F * exp(-1im * p.ωd * t)
    p = (F = F, ωd = ωd)

    # Time-dependent Hamiltonian as tuple (lab frame with drive)
    H_td = (H_lab, (a, coef1), (a', coef2))

    # Time-dependent Hamiltonian as QobjEvo
    H_td2 = QobjEvo(H_td)

    # Time-dependent Liouvillian
    L_td = liouvillian(H_td2)

    tlist_td = range(0, 10 / γ, 100)

    ## sesolve (time-dependent) ##

    SUITE["Time Evolution"]["time-dependent"]["sesolve"]["Tuple"] =
        @benchmarkable sesolve($H_td, $ψ0, $tlist_td, e_ops = $e_ops, progress_bar = Val(false), params = $p)

    SUITE["Time Evolution"]["time-dependent"]["sesolve"]["QobjEvo"] =
        @benchmarkable sesolve($H_td2, $ψ0, $tlist_td, e_ops = $e_ops, progress_bar = Val(false), params = $p)

    ## mesolve (time-dependent) ##

    SUITE["Time Evolution"]["time-dependent"]["mesolve"]["Tuple"] =
        @benchmarkable mesolve($H_td, $ψ0, $tlist_td, $c_ops, e_ops = $e_ops, progress_bar = Val(false), params = $p)

    SUITE["Time Evolution"]["time-dependent"]["mesolve"]["QobjEvo"] =
        @benchmarkable mesolve($H_td2, $ψ0, $tlist_td, $c_ops, e_ops = $e_ops, progress_bar = Val(false), params = $p)

    SUITE["Time Evolution"]["time-dependent"]["mesolve"]["Liouvillian"] =
        @benchmarkable mesolve($L_td, $ψ0, $tlist_td, $c_ops, e_ops = $e_ops, progress_bar = Val(false), params = $p)

    ## mcsolve (time-dependent) ##

    SUITE["Time Evolution"]["time-dependent"]["mcsolve"]["Tuple"]["Serial"] = @benchmarkable mcsolve(
        $H_td,
        $ψ0,
        $tlist_td,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        params = $p,
        ensemblealg = EnsembleSerial(),
    )

    SUITE["Time Evolution"]["time-dependent"]["mcsolve"]["Tuple"]["Multithreaded"] = @benchmarkable mcsolve(
        $H_td,
        $ψ0,
        $tlist_td,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        params = $p,
        ensemblealg = EnsembleThreads(),
    )

    SUITE["Time Evolution"]["time-dependent"]["mcsolve"]["QobjEvo"]["Serial"] = @benchmarkable mcsolve(
        $H_td2,
        $ψ0,
        $tlist_td,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        params = $p,
        ensemblealg = EnsembleSerial(),
    )

    SUITE["Time Evolution"]["time-dependent"]["mcsolve"]["QobjEvo"]["Multithreaded"] = @benchmarkable mcsolve(
        $H_td2,
        $ψ0,
        $tlist_td,
        $c_ops,
        ntraj = 100,
        e_ops = $e_ops,
        progress_bar = Val(false),
        params = $p,
        ensemblealg = EnsembleThreads(),
    )

    return nothing
end
