@testset "Block Diagonal Form" begin
    H, c_ops, a = driven_dissipative_kerr()
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
