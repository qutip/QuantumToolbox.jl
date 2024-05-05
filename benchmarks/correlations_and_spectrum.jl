function benchmark_correlations_and_spectrum()
    ωc = 1
    ωq = 1
    g = 0.05
    ωd = 0.95
    F = 0.1
    nth = 7
    γ = 0.1
    
    Δc = ωc - ωd
    Δq = ωq - ωd
    
    # Operators definition
    N = 50 # cutoff for the cavity Hilbert space
    a = tensor(destroy(N), qeye(2))
    σm = tensor(qeye(N), sigmam())
    σz = tensor(qeye(N), sigmaz())
    
    # Hamiltonian
    H = Δc * a' * a + Δq * σz / 2 + g * (a' * σm + a * σm') + F * (a + a')
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a', sqrt(γ) * σm]

    ω_l = range(0, 3, length=1000)

    SUITE["Correlations and Spectrum"]["FFT Correlation"] = @benchmarkable spectrum($H, $ω_l, $(a'), $a, $c_ops, solver=FFTCorrelation(), progress_bar=false)

    SUITE["Correlations and Spectrum"]["Exponential Series"] = @benchmarkable spectrum($H, $ω_l, $(a'), $a, $c_ops)
end

benchmark_correlations_and_spectrum()