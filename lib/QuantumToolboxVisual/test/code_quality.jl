@testset "Code quality (QuantumToolboxVisual)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxVisual; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxVisual; target_modules = (QuantumToolboxVisual,), ignore_missing_comparison = true)
    end
end
