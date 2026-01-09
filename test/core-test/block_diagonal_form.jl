@testitem "Block Diagonal Form" begin
    # Block Diagonal Form
    N = 20
    Δ = 0
    G = 5
    tg = 0
    θ = atan(tg)
    U = sin(θ)
    κ2 = cos(θ)
    κϕ = 1.0e-3
    nth = 0.0

    a = destroy(N)
    ad = create(N)
    H = -Δ * ad * a + G / 2 * (ad^2 + a^2) + U / 2 * (ad^2 * a^2)
    c_ops = [√(κ2) * a^2, √(κϕ) * ad * a]
    L = liouvillian(H, c_ops)

    bdf = block_diagonal_form(L)
    L_bd = bdf.B
    block_sizes = bdf.block_sizes
    blocks = bdf.blocks

    @test size(L_bd) == size(L)
    @test length(block_sizes) == 4
    @test length(blocks) == 4
    @test sum(block_sizes .== 100) == 4

    @testset "Type Inference (block_diagonal_form)" begin
        @inferred block_diagonal_form(L)
    end
end
