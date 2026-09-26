setprecision(128) # Instead of 256. This speeds up the tests.

using QuantumToolbox

include("setup.jl") # module TESetup (parameters and operators) are defined in this file

@testset "Arbitrary Precision (mesolve)" begin
    sol = mesolve(TESetup.H, TESetup.ψ0, TESetup.tlist, TESetup.c_ops; progress_bar = Val(false))
    sol_big = mesolve(TESetup.H_big, TESetup.ψ0_big, TESetup.tlist, TESetup.c_ops_big; progress_bar = Val(false))

    @test eltype(sol_big.states[1]) == Complex{BigFloat}

    # Test all fidelities are close to 1
    @test all((x) -> isapprox(hilbert_dist(x[1], x[2]), 0; atol = 1.0e-10), zip(sol.states, sol_big.states))

    @inferred mesolve(TESetup.H_big, TESetup.ψ0_big, TESetup.tlist, TESetup.c_ops_big; progress_bar = Val(false))
end
