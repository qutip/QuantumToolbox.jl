@testitem "CODATA Physical Constants" begin
    c = PhysicalConstants.c
    h = PhysicalConstants.h
    ħ = PhysicalConstants.ħ
    μ0 = PhysicalConstants.μ0
    ϵ0 = PhysicalConstants.ϵ0

    @test h / ħ ≈ 2 * π
    @test μ0 / (4.0e-7 * π) ≈ 1.0
    @test c^2 * μ0 * ϵ0 ≈ 1.0

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
        v = rand(Float64)
        _unit_list = [:J, :eV, :meV, :GHz, :mK]
        for u1 in _unit_list
            for u2 in _unit_list
                @inferred convert_unit(v, u1, u2)
            end
        end
    end
end
