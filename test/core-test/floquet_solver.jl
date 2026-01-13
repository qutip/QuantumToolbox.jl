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
    tlist = range(0.0, 2T, length=101)
    floquet_basis = FloquetBasis(H_tuple, T, tlist)

    
end