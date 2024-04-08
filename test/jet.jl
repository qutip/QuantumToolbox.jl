@testset "Code quality (JET.jl)" begin
    if VERSION >= v"1.9"
        JET.test_package(QuantumToolbox; target_defined_modules=true, ignore_missing_comparison=true)
    end
end