@testset "Permutation" begin
    # Block Diagonal Form
    N = 20
    Δ = 0
    G = 5
    tg = 0
    θ = atan(tg)
    U = sin(θ)
    κ2 = cos(θ)
    κ1 = 0.0
    κϕ = 1e-3
    nth = 0.0

    a = destroy(N)
    ad = create(N)
    H = -Δ * ad * a + G / 2 * (ad^2 + a^2) + U / 2 * (ad^2 * a^2)
    c_ops = [√(κ2) * a^2, √(κ1 * (nth + 1)) * a, √(κ1 * nth) * ad, √(κϕ) * ad * a]
    L = liouvillian(H, c_ops)

    P, L_bd, block_sizes = bdf(L)
    blocks_list, block_indices = get_bdf_blocks(L_bd, block_sizes)
    @test size(L_bd) == size(L)
    @test length(block_sizes) == 4
    @test length(blocks_list) == 4
    @test length(block_indices) == 4
    @test sum(block_sizes .== 100) == 4
end
