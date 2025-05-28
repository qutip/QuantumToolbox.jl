@testset "Makie Extension" verbose = true begin
    ψ = normalize(coherent(50, 5.0) + coherent(50, -5.0))
    xvec = yvec = -15.0:0.1:15.0
    wig = transpose(wigner(ψ, xvec, yvec))

    @test_throws ArgumentError plot_wigner(ψ; library = :Makie, xvec = xvec, yvec = yvec)

    @test_throws ArgumentError plot_fock_distribution(ψ; library = :Makie)

    using Makie

    fig, ax, hm =
        plot_wigner(ψ; library = Val(:Makie), xvec = xvec, yvec = yvec, projection = Val(:two_dim), colorbar = true)
    @test fig isa Figure
    @test ax isa Axis
    @test hm isa Heatmap
    @test all(isapprox.(hm[3].val, wig, atol = 1e-6))

    fig, ax, surf =
        plot_wigner(ψ; library = Val(:Makie), xvec = xvec, yvec = yvec, projection = Val(:three_dim), colorbar = true)
    @test fig isa Figure
    @test ax isa Axis3
    @test surf isa Surface
    @test all(isapprox.(surf[3].val, wig, atol = 1e-6))

    fig = Figure()
    pos = fig[2, 3]
    fig1, ax, hm = plot_wigner(
        ψ;
        library = Val(:Makie),
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
        library = Val(:Makie),
        xvec = xvec,
        yvec = yvec,
        projection = Val(:three_dim),
        colorbar = true,
        location = pos,
    )
    @test fig1 === fig
    @test fig[2, 3].layout.content[1].content[1, 1].layout.content[1].content === ax

    fig = Figure()
    pos = fig[2, 3]
    fig1, ax = plot_fock_distribution(ψ; library = Val(:Makie), location = pos)
    @test fig1 === fig
    @test fig[2, 3].layout.content[1].content[1, 1].layout.content[1].content === ax

    fig = Figure()
    pos = fig[2, 3]
    fig1, ax = @test_logs (:warn,) plot_fock_distribution(ψ * 2; library = Val(:Makie), location = pos)

    ρ = 0.7*ket2dm(basis(2, 0)) + 0.3*ket2dm(basis(2, 1))
    fig, ax = plot_bloch(ρ)
    @test fig isa Figure
    @test ax isa Axis3
    ψ = (basis(2, 0) + basis(2, 1))/√2
    fig, ax = plot_bloch(ψ)
    @test fig isa Figure
    @test ax isa Axis3
    ϕ = dag(ψ)
    fig, ax = plot_bloch(ϕ)
    @test fig isa Figure
    @test ax isa Axis3
    fig = Figure()
    pos = fig[1, 1]
    fig1, ax = plot_bloch(ψ; location = pos)
    @test fig1 === fig

    b = Bloch()
    add_points!(b, [0.0, 0.0, 1.0])
    @test length(b.points) == 1
    @test b.points[1] == [0.0, 0.0, 1.0]
    pts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    add_points!(b, pts)
    @test length(b.points) == 3
    @test b.points[2] == pts[1]
    b = Bloch()
    add_points!(b, [0.0, 0.0, 1.0])
    @test length(b.points) == 1
    @test b.points[1] == [0.0, 0.0, 1.0]
    pts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    add_points!(b, pts)
    @test length(b.points) == 3
    @test b.points[2] == pts[1]
    b = Bloch()
    add_vectors!(b, [1.0, 1.0, 0.0])
    @test length(b.vectors) == 1
    @test isapprox(norm(b.vectors[1]), 1.0)
    vecs = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
    add_vectors!(b, vecs)
    @test length(b.vectors) == 3
    @test all(norm(v) ≈ 1.0 for v in b.vectors)
    b = Bloch()
    add_line!(b, [0, 0, 0], [1, 0, 0])
    @test length(b.lines) == 1
    @test b.lines[1][2] == [1.0, 0.0, 0.0]
    add_arc!(b, [0, 0, 1], [0, 1, 0], [1, 0, 0])
    @test length(b.arcs) == 1
    @test b.arcs[1][3] == [1.0, 0.0, 0.0]
    b = Bloch()
    add_points!(b, [0.0, 0.0, 1.0])
    add_vectors!(b, [1.0, 0.0, 0.0])
    add_line!(b, [0, 0, 0], [1, 0, 0])
    add_arc!(b, [0, 1, 0], [0, 0, 1], [1, 0, 0])
    clear!(b)
    @test isempty(b.points)
    @test isempty(b.vectors)
    @test isempty(b.lines)
    @test isempty(b.arcs)
    b = Bloch()
    add_points!(b, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    add_vectors!(b, [[1, 1, 0], [0, 1, 1]])
    add_line!(b, [0, 0, 0], [1, 1, 1])
    add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
    try
        fig, ax = QuantumToolbox.render(b)
        @test !isnothing(fig)
        @test !isnothing(ax)
    catch e
        @test false
        @info "Render threw unexpected error" exception=e
    end
end
