function benchmark_correlations_and_spectrum()
    N = 15
    ω = 1
    γ = 0.1
    nth = 0.02

    a = destroy(N)
    H = ω * a' * a
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a']

    ω_l = range(0, 3, length = 1000)

    SUITE["Correlations and Spectrum"]["FFT Correlation"] =
        @benchmarkable spectrum($H, $ω_l, $(a'), $a, $c_ops, solver = FFTCorrelation(), progress_bar = false)

    return SUITE["Correlations and Spectrum"]["Exponential Series"] =
        @benchmarkable spectrum($H, $ω_l, $(a'), $a, $c_ops)
end

benchmark_correlations_and_spectrum()
