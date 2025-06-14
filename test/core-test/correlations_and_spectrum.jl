@testitem "Correlations and Spectrum" begin
    N = 10
    Id = qeye(N)
    a = destroy(N)
    H = a' * a
    c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']

    t_l = range(0, 333 * π, length = 1000)
    corr1 = correlation_2op_1t(H, nothing, t_l, c_ops, a', a; progress_bar = Val(false))
    corr2 = correlation_3op_1t(H, nothing, t_l, c_ops, Id, a', a; progress_bar = Val(false))
    ω_l1, spec1 = spectrum_correlation_fft(t_l, corr1)

    ω_l2 = range(0, 3, length = 1000)
    spec2 = spectrum(H, ω_l2, c_ops, a', a)
    spec3 = spectrum(H, ω_l2, c_ops, a', a; solver = PseudoInverse())
    spec4 = spectrum(H, ω_l2, c_ops, a', a; solver = Lanczos())

    spec1 = spec1 ./ maximum(spec1)
    spec2 = spec2 ./ maximum(spec2)
    spec3 = spec3 ./ maximum(spec3)
    spec4 = spec4 ./ maximum(spec4)

    test_func1 = maximum(real.(spec1)) * (0.1 / 2)^2 ./ ((ω_l1 .- 1) .^ 2 .+ (0.1 / 2)^2)
    test_func2 = maximum(real.(spec2)) * (0.1 / 2)^2 ./ ((ω_l2 .- 1) .^ 2 .+ (0.1 / 2)^2)
    idxs1 = test_func1 .> 0.05
    idxs2 = test_func2 .> 0.05
    @test sum(abs2.(spec1[idxs1] .- test_func1[idxs1])) / sum(abs2.(test_func1[idxs1])) < 0.01
    @test sum(abs2.(spec2[idxs2] .- test_func2[idxs2])) / sum(abs2.(test_func2[idxs2])) < 0.01
    @test all(corr1 .≈ corr2)
    @test all(spec2 .≈ spec3)
    @test all(spec2 .≈ spec4)

    @testset "Type Inference spectrum" begin
        @inferred correlation_2op_1t(H, nothing, t_l, c_ops, a', a; progress_bar = Val(false))
        @inferred spectrum_correlation_fft(t_l, corr1)
        @inferred spectrum(H, ω_l2, c_ops, a', a)
        @inferred spectrum(H, ω_l2, c_ops, a', a; solver = PseudoInverse())
        @inferred spectrum(H, ω_l2, c_ops, a', a; solver = Lanczos())
    end

    @testset "Verbose mode Lanczos" begin
        cout = stdout
        r, w = redirect_stdout()
        nout = @async read(r, String)
        spectrum(H, ω_l2, c_ops, a', a; solver = Lanczos(verbose = 2, maxiter = 2, tol = 1e-16));
        redirect_stdout(cout)
        close(w)
        out = fetch(nout)
        outlines = split(out, '\n', keepempty = false)
        @test last(outlines) == "spectrum(): Consider increasing maxiter and/or tol"
    end

    # tlist and τlist checks
    t_fft_wrong = [0, 1, 10]
    t_wrong1 = [1, 2, 3]
    t_wrong2 = [-1, 0, 1]
    @test_throws ArgumentError spectrum_correlation_fft(t_fft_wrong, corr1)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_l, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_l, t_wrong2, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_l, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_l, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_wrong2, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_wrong2, c_ops, Id, a', a)

    @testset "Deprecated Errors" begin
        ρ0 = rand_dm(N)
        @test_throws ErrorException FFTCorrelation()
        @test_throws ErrorException correlation_3op_2t(H, ρ0, t_l, t_l, a, a', a, c_ops)
        @test_throws ErrorException correlation_3op_1t(H, ρ0, t_l, a, a', a, c_ops)
        @test_throws ErrorException correlation_2op_2t(H, ρ0, t_l, t_l, a', a, c_ops)
        @test_throws ErrorException correlation_2op_1t(H, ρ0, t_l, a', a, c_ops)
    end
end
