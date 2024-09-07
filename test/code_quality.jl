using Aqua, JET

@testset "Code quality" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolbox; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolbox; target_defined_modules = true, ignore_missing_comparison = true)
    end
end
