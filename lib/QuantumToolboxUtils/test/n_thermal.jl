@testitem "n_thermal" begin
    ω1 = rand(Float64)
    ω2 = rand(Float64)
    @test n_thermal(0, ω2) == 0.0
    @test n_thermal(ω1, 0) == 0.0
    @test n_thermal(ω1, -ω2) == 0.0
    @test n_thermal(ω1, ω2) == 1 / (exp(ω1 / ω2) - 1)
    @test typeof(n_thermal(Int32(2), Int32(3))) == Float32
    @test typeof(n_thermal(Float32(2), Float32(3))) == Float32
    @test typeof(n_thermal(Int64(2), Int32(3))) == Float64
    @test typeof(n_thermal(Int32(2), Int64(3))) == Float64
    @test typeof(n_thermal(Float64(2), Float32(3))) == Float64
    @test typeof(n_thermal(Float32(2), Float64(3))) == Float64

    @testset "Type Inference" begin
        v = rand(Float64)
        @inferred n_thermal(v, Int32(123))
    end
end
