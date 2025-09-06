@testset "Makie Extension" verbose = true begin
    ψ = normalize(coherent(50, 5.0) + coherent(50, -5.0))
    xvec = yvec = -15.0:0.1:15.0
    wig = transpose(wigner(ψ, xvec, yvec))

    # Makie unload errors
    @test_throws ArgumentError plot_wigner(ψ; library = :Makie, xvec = xvec, yvec = yvec)
    @test_throws ArgumentError plot_fock_distribution(ψ; library = :Makie)
    @test_throws ArgumentError plot_bloch(ψ; library = :Makie)

    using Makie

    fig, ax, hm =
        plot_wigner(ψ; library = Val(:Makie), xvec = xvec, yvec = yvec, projection = Val(:two_dim), colorbar = true)
    @test fig isa Figure
    @test ax isa Axis
    @test hm isa Heatmap
    @test all(isapprox.(hm[3].value.x, wig, atol = 1e-6))

    fig, ax, surf =
        plot_wigner(ψ; library = Val(:Makie), xvec = xvec, yvec = yvec, projection = Val(:three_dim), colorbar = true)
    @test fig isa Figure
    @test ax isa Axis3
    @test surf isa Surface
    @test all(isapprox.(surf[3].value.x, wig, atol = 1e-6))

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
end

@testset "Makie Bloch sphere" begin
    ρ = 0.7 * ket2dm(basis(2, 0)) + 0.3 * ket2dm(basis(2, 1))
    fig, lscene = plot_bloch(ρ)
    @test fig isa Figure
    @test lscene isa LScene

    ψ = (basis(2, 0) + basis(2, 1)) / √2
    fig, lscene = plot_bloch(ψ)
    @test fig isa Figure
    @test lscene isa LScene

    ϕ = dag(ψ)
    fig, lscene = plot_bloch(ϕ)
    @test fig isa Figure
    @test lscene isa LScene

    fig = Figure()
    ax = Axis(fig[1, 1])
    pos = fig[1, 2]
    fig1, lscene = plot_bloch(ψ; location = pos)
    @test fig1 === fig

    b = Bloch()
    add_points!(b, [0.0, 0.0, 1.0])
    @test length(b.points) == 1
    @test b.points[1] ≈ [0.0, 0.0, 1.0]

    pts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    add_points!(b, hcat(pts...))
    @test length(b.points) == 2
    @test b.points[2] ≈ hcat(pts...)
    @test_throws ArgumentError add_points!(b, [1 2 3 4])
    @test_throws ArgumentError add_points!(b, pts; meth = :wrong)

    b = Bloch()
    add_vectors!(b, [1.0, 1.0, 0.0])
    @test length(b.vectors) == 1
    @test isapprox(norm(b.vectors[1]), √2)

    vecs = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
    add_vectors!(b, vecs)
    @test length(b.vectors) == 3
    @test isapprox(norm(b.vectors[2]), 1.0)
    @test isapprox(norm(b.vectors[3]), 1.0)

    vec_correct = [1, 0, 0]
    vec_wrong = [1, 0]
    b = Bloch()
    add_line!(b, [0, 0, 0], [1, 0, 0])
    @test length(b.lines) == 1
    @test b.lines[1][1][3] ≈ [0.0, 0.0]
    @test_throws ArgumentError add_line!(b, vec_wrong, vec_correct)
    @test_throws ArgumentError add_line!(b, vec_correct, vec_wrong)

    add_arc!(b, [0, 0, 1], [0, 1, 0], [1, 0, 0])
    @test length(b.arcs) == 1
    @test b.arcs[1][3] == [1.0, 0.0, 0.0]
    @test_throws ArgumentError add_arc!(b, vec_wrong, vec_correct)
    @test_throws ArgumentError add_arc!(b, vec_correct, vec_wrong)
    @test_throws ArgumentError add_arc!(b, vec_wrong, vec_correct, vec_correct)
    @test_throws ArgumentError add_arc!(b, vec_correct, vec_wrong, vec_correct)
    @test_throws ArgumentError add_arc!(b, vec_correct, vec_correct, vec_wrong)

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
    add_points!(b, hcat([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    add_vectors!(b, [[1, 1, 0], [0, 1, 1]])
    add_line!(b, [0, 0, 0], [1, 1, 1])
    add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
    try
        fig, lscene = render(b)
        @test !isnothing(fig)
        @test !isnothing(lscene)
    catch e
        @test false
        @info "Render threw unexpected error" exception=e
    end

    b = Bloch()
    ψ₁ = normalize(basis(2, 0) + basis(2, 1))
    ψ₂ = normalize(basis(2, 0) - im * basis(2, 1))
    add_line!(b, ψ₁, ψ₂; fmt = "r--")
    try
        fig, lscene = render(b)
        @test !isnothing(fig)
        @test !isnothing(lscene)
    catch e
        @test false
        @info "Render threw unexpected error" exception=e
    end

    # test `state to Bloch vector` conversion and `add_states!` function
    b = Bloch()
    Pauli_Ops = [sigmax(), sigmay(), sigmaz()]
    ψ = rand_ket(2)
    ρ = rand_dm(2)
    states = [ψ, ρ]
    x = basis(2, 0) + basis(2, 1)             # unnormalized Ket
    ρ1 = 0.3 * rand_dm(2) + 0.4 * rand_dm(2)  # unnormalized density operator
    ρ2 = Qobj(rand(ComplexF64, 2, 2))         # unnormalized and non-Hermitian Operator
    add_states!(b, states, kind = :vector)
    add_states!(b, states, kind = :point)
    @test length(b.vectors) == 2
    @test length(b.points) == 1
    @test all(expect(Pauli_Ops, ψ) .≈ (b.vectors[1]))
    @test all(expect(Pauli_Ops, ρ) .≈ (b.vectors[2]))
    @test all([b.vectors[j][k] ≈ b.points[1][k, j] for j in (1, 2) for k in (1, 2, 3)])
    @test_logs (:warn,) (:warn,) (:warn,) (:warn,) add_states!(b, [x, ρ1, ρ2])
    @test length(b.vectors) == 5
    @test_throws ArgumentError add_states!(b, states, kind = :wrong)

    th = range(0, 2π; length = 20)
    xp = cos.(th);
    yp = sin.(th);
    zp = zeros(20);
    pnts = [xp, yp, zp];
    pnts = Matrix(hcat(xp, yp, zp)');
    add_points!(b, pnts);
    vec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    add_vectors!(b, vec);
    add_line!(b, [1, 0, 0], [0, 1, 0])
    add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
    try
        fig, lscene = render(b)
        @test !isnothing(fig)
        @test !isnothing(lscene)
    catch e
        @test false
        @info "Render threw unexpected error" exception=e
    end

    b = Bloch()
    ψ₁ = normalize(basis(2, 0) + basis(2, 1))
    ψ₂ = normalize(basis(2, 0) - im * basis(2, 1))
    ψ₃ = basis(2, 0)
    add_line!(b, ψ₁, ψ₂; fmt = "r--")
    add_arc!(b, ψ₁, ψ₂)
    add_arc!(b, ψ₂, ψ₃, ψ₁)
    add_states!(b, [ψ₂, ψ₃], kind = :point, meth = :l)
    try
        fig, lscene = render(b)
        @test !isnothing(fig)
        @test !isnothing(lscene)
    catch e
        @test false
        @info "Render threw unexpected error" exception=e
    end

    # if render location is given as lscene, should return the same Figure and LScene
    b = Bloch()
    fig1, lscene1 = render(b)
    add_states!(b, ψ₁)
    fig2, lscene2 = render(b, location = lscene1)
    @test fig2 === fig1
    @test lscene2 === lscene1
end
