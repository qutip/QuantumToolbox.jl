# Setup problem for arbitrary precision tests
N = 20
Δ = 1.0
U = 0.1
F = 0.1
κ = 0.1

a = destroy(N)
ψ0 = fock(N, 1)
H = Δ * a' * a + U / 2 * a' * a' * a * a + F * (a + a')

c_ops = [sqrt(κ) * a]

tlist = range(0, 40, 100)

a_big = destroy(Complex{BigFloat}, N)
ψ0_big = fock(Complex{BigFloat}, N, 1)
H_big = Δ * a_big' * a_big + U / 2 * a_big' * a_big' * a_big * a_big + F * (a_big + a_big')

c_ops_big = [sqrt(BigFloat(κ)) * a_big]
