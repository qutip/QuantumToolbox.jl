@testset "Entanglement" begin
    g = fock(2, 1)
    e = fock(2, 0)
    state = normalize(kron(g, e) + kron(e, g))
    rho = state * state'
    @test entanglement(state, 1) / log(2) ≈ 1
    @test entanglement(rho, 1) / log(2) ≈ 1
end
