using QuantumToolbox
using ForwardDiff

@testset "ForwardDiff for thermal_dm" begin
    N = 100  # dimension of the system
    N_op = num(N) # number operator

    function N_expect(p) # p = [n]
        ρT = thermal_dm(N, p[1]; sparse = Val(true))
        return real(expect(N_op, ρT))
    end

    # Average photon number for thermal state
    n = rand(Float64)
    p = [n]

    # Use ForwardDiff.gradient to compute the gradient
    grad = ForwardDiff.gradient(N_expect, p)[1]

    # Compare the result
    @test isapprox(grad, 1.0; atol = 1.0e-6)
end