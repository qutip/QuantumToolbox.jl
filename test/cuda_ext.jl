using CUDA

CUDA.versioninfo()

@testset "CUDA Extension" begin
    @test_throws DomainError cu(sigmax(); word_size=16)
end