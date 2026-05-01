@testset "Code quality (QuantumToolboxUtils)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxUtils; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxUtils; target_modules = (QuantumToolboxUtils,), ignore_missing_comparison = true)
    end
end
