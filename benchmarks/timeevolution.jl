function benchmark_timeevolution()
    ωc = 1
    ωq = 1
    g = 0.05
    ωd = 0.95
    F = 0.1
    
    Δc = ωc - ωd
    Δq = ωq - ωd
    
    # Operators definition
    N = 50 # cutoff for the cavity Hilbert space
    a = tensor(destroy(N), qeye(2))
    σm = tensor(qeye(N), sigmam())
    σz = tensor(qeye(N), sigmaz())
    
    # Hamiltonian
    H = Δc * a' * a + Δq * σz / 2 + g * (a' * σm + a * σm') + F * (a + a')
    
    e_ops = [a'*a, σz]
    
    # Initial state
    ψ0 = tensor(coherent(N, 0), fock(2, 1))

    ## sesolve ##
    
    tlist = range(0, 2π * 10 / g, 1000)
    
    SUITE["Time Evolution"]["time-independent"]["sesolve"] = @benchmarkable sesolve($H, $ψ0, $tlist, e_ops=$e_ops, progress_bar=false)

    ## mesolve ##

    nth = 7
    γ = 0.1
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a', sqrt(γ) * σm]
    
    tlist = range(0, 10/γ, 100)
    
    SUITE["Time Evolution"]["time-independent"]["mesolve"] = @benchmarkable mesolve($H, $ψ0, $tlist, $c_ops, e_ops=$e_ops, progress_bar=false)

    ## mcsolve ##

    ntraj = 100

    SUITE["Time Evolution"]["time-independent"]["mcsolve"] = @benchmarkable mcsolve($H, $ψ0, $tlist, $c_ops, n_traj=ntraj, e_ops=$e_ops)
end

benchmark_timeevolution()