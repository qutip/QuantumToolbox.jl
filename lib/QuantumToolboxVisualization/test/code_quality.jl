@testset "Code quality (QuantumToolboxVisualization)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxVisualization; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxVisualization; target_modules = (QuantumToolboxVisualization,), ignore_missing_comparison = true)
    end
end
