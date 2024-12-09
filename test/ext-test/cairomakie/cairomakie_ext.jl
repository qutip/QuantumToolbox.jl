@testset "CairoMakie Extension" verbose = true begin
    ψ = normalize(coherent(50, 5.0) + coherent(50, -5.0))
    xvec = yvec = -15.0:0.1:15.0
    wig = wigner(ψ, xvec, yvec)'
    
    @test_throws ArgumentError plot_wigner(ψ; library = :CairoMakie, xvec = xvec, yvec = yvec)

    using CairoMakie

    fig, ax, hm = plot_wigner(
        ψ;
        library = Val(:CairoMakie),
        xvec = xvec,
        yvec = yvec,
        projection = Val(:two_dim),
        colorbar = true,
    )
    @test fig isa Figure
    @test ax isa Axis
    @test hm isa Heatmap
    @test all(isapprox.(hm[3].val, wig, atol = 1e-6))

    fig, ax, surf = plot_wigner(
        ψ;
        library = Val(:CairoMakie),
        xvec = xvec,
        yvec = yvec,
        projection = Val(:three_dim),
        colorbar = true,
    )
    @test fig isa Figure
    @test ax isa Axis3
    @test surf isa Surface
    @test all(isapprox.(surf[3].val, wig, atol = 1e-6))

    fig = Figure()
    pos = fig[2, 3]
    fig1, ax, hm = plot_wigner(
        ψ;
        library = Val(:CairoMakie),
        xvec = xvec,
        yvec = yvec,
        projection = Val(:two_dim),
        colorbar = true,
        location = pos,
    )
    @test fig1 === fig
    @test fig[2, 3].layout.content[1].content[1, 1].layout.content[1].content === ax

    fig = Figure()
    pos = fig[2, 3]
    fig1, ax, surf = plot_wigner(
        ψ;
        library = Val(:CairoMakie),
        xvec = xvec,
        yvec = yvec,
        projection = Val(:three_dim),
        colorbar = true,
        location = pos,
    )
    @test fig1 === fig
    @test fig[2, 3].layout.content[1].content[1, 1].layout.content[1].content === ax
end
