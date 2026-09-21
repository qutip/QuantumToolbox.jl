setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox

include("setup.jl") # parameters and operators are defined in this file

@testset "Arbitrary Precision (mesolve)" begin
    sol = mesolve(H, ψ0, tlist, c_ops; progress_bar = Val(false))
    sol_big = mesolve(H_big, ψ0_big, tlist, c_ops_big; progress_bar = Val(false))

    @test eltype(sol_big.states[1]) == Complex{BigFloat}

    # Test all fidelities are close to 1
    @test all((x) -> isapprox(hilbert_dist(x[1], x[2]), 0; atol = 1.0e-10), zip(sol.states, sol_big.states))

    @inferred mesolve(H_big, ψ0_big, tlist, c_ops_big; progress_bar = Val(false))
end
