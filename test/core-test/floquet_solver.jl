using QuantumToolbox
using Test 
using Random

@testitem "Floquet Solver" begin
    N = 10
    a = destroy(N)
    a_d = a'
    coef(p,t) = cos(t)
    H0 = num(N)
    H1 = a + a_d
    H_tuple = (H0,(H1,coef))
    T = 2 * π
    tlist = range(0.0, 3T, length=101)
    floquet_basis = FloquetBasis(H_tuple, T, tlist)
    psi0 = rand_ket(N)
    floquet_psi0 = to_floquet_basis(floquet_basis, psi0)
    sol = sesolve(H_tuple, psi0, tlist, e_ops = [], saveat = tlist)
    states = sol.states
    for (t,state) in zip(tlist,states)
        from_floquet = from_floquet_basis(floquet_basis, floquet_psi0, t)
        @test overlap(state, from_floquet) ≈ 1.0 atol=8e-5  
end