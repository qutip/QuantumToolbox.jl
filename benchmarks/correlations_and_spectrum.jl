function _calculate_fft_spectrum(H, tlist, c_ops, A, B)
    corr = correlation_2op_1t(H, nothing, tlist, c_ops, A, B; progress_bar = Val(false))
    ωlist, spec = spectrum_correlation_fft(tlist, corr)
    return nothing
end

function benchmark_correlations_and_spectrum!(SUITE)
    N = 15
    ω = 1
    γ = 0.1
    nth = 0.02

    a = destroy(N)
    H = ω * a' * a
    c_ops = [sqrt(γ * (nth + 1)) * a, sqrt(γ * nth) * a']

    ω_l = range(0, 3, length = 1000)
    t_l = range(0, 333 * π, length = 1000)

    PI_solver = PseudoInverse()

    L_solver = Lanczos()

    SUITE["Correlations and Spectrum"]["FFT Correlation"] =
        @benchmarkable _calculate_fft_spectrum($H, $t_l, $c_ops, $(a'), $a)

    SUITE["Correlations and Spectrum"]["Spectrum"]["Exponential Series"] =
        @benchmarkable spectrum($H, $ω_l, $c_ops, $(a'), $a)

    SUITE["Correlations and Spectrum"]["Spectrum"]["Pseudo Inverse"] =
        @benchmarkable spectrum($H, $ω_l, $c_ops, $(a'), $a, solver = $PI_solver)

    SUITE["Correlations and Spectrum"]["Spectrum"]["Lanczos"] =
        @benchmarkable spectrum($H, $ω_l, $c_ops, $(a'), $a, solver = $L_solver)

    return nothing
end
