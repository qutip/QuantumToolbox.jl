@testset "Utilities" verbose = true begin
    @testset "n_thermal" begin
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
    end

    @testset "CODATA Physical Constants" begin
        c = PhysicalConstants.c
        h = PhysicalConstants.h
        ħ = PhysicalConstants.ħ
        μ0 = PhysicalConstants.μ0
        ϵ0 = PhysicalConstants.ϵ0

        @test h / ħ ≈ 2 * π
        @test μ0 / (4e-7 * π) ≈ 1.0
        @test c^2 * μ0 * ϵ0 ≈ 1.0
    end

    @testset "convert unit" begin
        V = 100 * rand(Float64)
        _unit_list = [:J, :eV, :meV, :MHz, :GHz, :K, :mK]
        for origin in _unit_list
            for middle in _unit_list
                for target in _unit_list
                    V_middle = convert_unit(V, origin, middle)
                    V_target = convert_unit(V_middle, middle, target)
                    V_origin = convert_unit(V_target, target, origin)
                    @test V ≈ V_origin
                end
            end
        end
        @test_throws ArgumentError convert_unit(V, :bad_unit, :J)
        @test_throws ArgumentError convert_unit(V, :J, :bad_unit)
    end

    @testset "Type Inference" begin
        v1 = rand(Float64)
        @inferred n_thermal(v1, Int32(123))

        _unit_list = [:J, :eV, :meV, :GHz, :mK]
        for u1 in _unit_list
            for u2 in _unit_list
                @inferred convert_unit(v1, u1, u2)
            end
        end
    end
end
