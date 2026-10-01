using Test
using QuantumToolbox

# module TESetup (parameters and operators) are defined in this file
# note that this file also sets the precision to 128 bits
include("setup.jl")

@testset "Arbitrary Precision (sesolve)" begin
    sol = sesolve(TESetup.H, TESetup.ψ0, TESetup.tlist; progress_bar = Val(false))
    sol_big = sesolve(TESetup.H_big, TESetup.ψ0_big, TESetup.tlist; progress_bar = Val(false))

    @test eltype(sol_big.states[1]) == Complex{BigFloat}

    # Test all fidelities are close to 1
    @test all((x) -> isapprox(fidelity(x[1], x[2]), 1; atol = 1.0e-7), zip(sol.states, sol_big.states))

    @inferred sesolve(TESetup.H_big, TESetup.ψ0_big, TESetup.tlist; progress_bar = Val(false))
end
