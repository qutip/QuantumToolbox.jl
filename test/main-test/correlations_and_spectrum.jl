#!PTR_MULTITHREAD
using Test
using QuantumToolbox

@testset "Correlations and Spectrum" begin
    @test Base.Threads.nthreads() > 1 # make sure the test runs with multi-threading

    N = 10
    Id = qeye(N)
    a = destroy(N)
    H = a' * a
    c_ops = [sqrt(0.1 * (0.01 + 1)) * a, sqrt(0.1 * (0.01)) * a']
    ρss = steadystate(H, c_ops)

    # test correlation_3op_2t (calculate g² correlations)
    Nt = 10
    Nτ = 8
    tlist = range(0, Nt, length = Nt)
    τlist = range(0, Nτ, length = Nτ)
    time_cutoff = 5.0 # exclude the cases where t + τ > time_cutoff
    A = a'
    B = a' * a
    C = a
    g2_corr1 = correlation_3op_2t(H, ρss, tlist, τlist, c_ops, A, B, C; progress_bar = Val(false), matrix_form = Val(false), max_t_plus_tau = nothing)
    g2_corr2 = correlation_3op_2t(H, ρss, tlist, τlist, c_ops, A, B, C; progress_bar = Val(false), matrix_form = Val(true), max_t_plus_tau = nothing)
    g2_cut_1 = correlation_3op_2t(H, ρss, tlist, τlist, c_ops, A, B, C; progress_bar = Val(false), matrix_form = Val(false), max_t_plus_tau = time_cutoff)
    g2_cut_2 = correlation_3op_2t(H, ρss, tlist, τlist, c_ops, A, B, C; progress_bar = Val(false), matrix_form = Val(true), max_t_plus_tau = time_cutoff)
    @test size(g2_corr1) == size(g2_corr2) == size(g2_cut_1) == size(g2_cut_2) == (Nt, Nτ)
    @test all(g2_corr1 .≈ g2_corr2)
    @test all(g2_cut_1 .≈ g2_cut_2)
    for (i, t) in enumerate(tlist)
        for (j, τ) in enumerate(τlist)
            if t + τ > time_cutoff
                @test g2_cut_1[i, j] == g2_cut_2[i, j] == 0 # exact zero
            else
                # atol accounts for the exact t + τ == cutoff boundary, where the adaptive solver takes a
                # slightly different step sequence depending on whether it's the trajectory's own endpoint (cutoff boundary) or an interior point (full computation)
                @test isapprox(g2_cut_1[i, j], g2_corr1[i, j]; atol = 1.0e-5)
                @test isapprox(g2_cut_2[i, j], g2_corr1[i, j]; atol = 1.0e-5)
            end
        end
    end

    # consistency check between correlation function and spectrum
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

    if VERSION >= v"1.11" # eigen is type unstable on v1.10
        @testset "Type Inference spectrum" begin
            @inferred correlation_2op_1t(H, nothing, t_l, c_ops, a', a; progress_bar = Val(false))
            @inferred spectrum_correlation_fft(t_l, corr1)
            @inferred spectrum(H, ω_l2, c_ops, a', a)
            @inferred spectrum(H, ω_l2, c_ops, a', a; solver = PseudoInverse())
            @inferred spectrum(H, ω_l2, c_ops, a', a; solver = Lanczos())
        end
    end

    @testset "Verbose mode Lanczos" begin
        cout = stdout
        r, w = redirect_stdout()
        nout = @async read(r, String)
        spectrum(H, ω_l2, c_ops, a', a; solver = Lanczos(verbose = 2, maxiter = 2, tol = 1.0e-16))
        redirect_stdout(cout)
        close(w)
        out = fetch(nout)
        outlines = split(out, '\n', keepempty = false)
        @test last(outlines) == "spectrum(): Consider increasing maxiter and/or tol"
    end

    @testset "Orthogonal input vectors Lanczos" begin
        @test_throws AssertionError spectrum(H, ω_l2, [c_ops[1]], a', a; solver = Lanczos())
    end

    # tlist and τlist checks
    t_fft_wrong = [0, 1, 10]
    t_wrong1 = [1, 2, 3]
    t_wrong2 = [-1, 0, 1]
    max_t_plus_tau_wrong = -1.0
    H_t = QobjEvo(H, (p, t) -> t)
    @test_throws ArgumentError spectrum_correlation_fft(t_fft_wrong, corr1)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_l, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_l, t_wrong2, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_l, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_l, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong1, t_wrong2, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_wrong1, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_wrong2, t_wrong2, c_ops, Id, a', a)
    @test_throws ArgumentError correlation_3op_2t(H, nothing, t_l, t_l, c_ops, Id, a', a; matrix_form = Val(true)) # requires specifying initial state
    @test_throws ArgumentError correlation_3op_2t(H_t, nothing, t_l, t_l, c_ops, Id, a', a)  # requires specifying initial state

    @testset "Deprecated Errors and Warnings" begin
        ρ0 = rand_dm(N)
        t_l = [0.0, 1.0] # make time list shorter
        @test_throws ErrorException FFTCorrelation()
        @test_logs (:warn,) correlation_3op_2t(H, ρ0, t_l, t_l, a, a', a, c_ops; progress_bar = Val(false))
        @test_logs (:warn,) correlation_3op_1t(H, ρ0, t_l, a, a', a, c_ops; progress_bar = Val(false))
        @test_logs (:warn,) correlation_2op_2t(H, ρ0, t_l, t_l, a', a, c_ops; progress_bar = Val(false))
        @test_logs (:warn,) correlation_2op_1t(H, ρ0, t_l, a', a, c_ops; progress_bar = Val(false))
    end
end
