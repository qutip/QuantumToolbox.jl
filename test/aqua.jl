using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(QuantumToolbox; ambiguities = false)
end
