using JET

@testset "Code quality (JET.jl)" begin
    JET.test_package(QuantumToolbox; target_defined_modules = true, ignore_missing_comparison = true)
end
