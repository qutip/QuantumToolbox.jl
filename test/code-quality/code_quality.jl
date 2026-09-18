@testset "Code quality (QuantumToolboxCore)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxCore; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxCore; target_modules = (QuantumToolboxCore,), ignore_missing_comparison = true)
    end
end

@testset "Code quality (QuantumToolboxVisual)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolboxVisual; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolboxVisual; target_modules = (QuantumToolboxVisual,), ignore_missing_comparison = true)
    end
end

@testset "Code quality (QuantumToolbox)" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolbox; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumToolbox; target_modules = (QuantumToolbox,), ignore_missing_comparison = true)
    end
end
