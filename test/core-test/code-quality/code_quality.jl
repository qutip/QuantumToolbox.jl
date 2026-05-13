@testset "Code quality" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolbox; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolbox; target_modules = (QuantumToolbox,), ignore_missing_comparison = true)
    end
end
