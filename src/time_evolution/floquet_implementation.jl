# Replicate the plots in Qutip Floquet tutorial

# First box

δ = 0.2 * 2 * π
ϵ0 = 1.0 * 2 * π
A = 2.5 * 2 * π
ω = 1.0 * 2 * π
H0 = - δ/2.0 * sigmax() - ϵ0/2.0 * sigmaz()
H1 = A/2.0 * sigmaz()
H = (H0, (H1, t -> cos(ω * t)))

# Second box

T = 2π / ω
floquet_basis = FloquetBasis(H,T)
