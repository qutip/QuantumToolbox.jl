@testset "Correlations and Spectrum" begin
    a = destroy(10)
    H = a' * a
    c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']

    ω_l = range(0, 3, length = 1000)
    ω_l1, spec1 = spectrum(H, ω_l, a', a, c_ops, solver = FFTCorrelation(), progress_bar = false)
    ω_l2, spec2 = spectrum(H, ω_l, a', a, c_ops)
    spec1 = spec1 ./ maximum(spec1)
    spec2 = spec2 ./ maximum(spec2)

    test_func1 = maximum(real.(spec1)) * (0.1 / 2)^2 ./ ((ω_l1 .- 1) .^ 2 .+ (0.1 / 2)^2)
    test_func2 = maximum(real.(spec2)) * (0.1 / 2)^2 ./ ((ω_l2 .- 1) .^ 2 .+ (0.1 / 2)^2)
    idxs1 = test_func1 .> 0.05
    idxs2 = test_func2 .> 0.05
    @test sum(abs2.(spec1[idxs1] .- test_func1[idxs1])) / sum(abs2.(test_func1[idxs1])) < 0.01
    @test sum(abs2.(spec2[idxs2] .- test_func2[idxs2])) / sum(abs2.(test_func2[idxs2])) < 0.01
end
