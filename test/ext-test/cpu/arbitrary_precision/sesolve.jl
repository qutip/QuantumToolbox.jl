setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox

include("setup.jl") # parameters and operators are defined in this file

@testset "Arbitrary Precision (sesolve)" begin
    sol = sesolve(H, ψ0, tlist; progress_bar = Val(false))
    sol_big = sesolve(H_big, ψ0_big, tlist; progress_bar = Val(false))

    @test eltype(sol_big.states[1]) == Complex{BigFloat}

    # Test all fidelities are close to 1
    @test all((x) -> isapprox(fidelity(x[1], x[2]), 1; atol = 1.0e-7), zip(sol.states, sol_big.states))

    @inferred sesolve(H_big, ψ0_big, tlist; progress_bar = Val(false))
end
