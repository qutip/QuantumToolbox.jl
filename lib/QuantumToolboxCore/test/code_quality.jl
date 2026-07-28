@testset "Code quality (QuantumToolboxCore)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxCore; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxCore; target_modules = (QuantumToolboxCore,), ignore_missing_comparison = true)
    end
end
