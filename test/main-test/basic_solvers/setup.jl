# Setup problem for arbitrary precision tests
module TESetup
using QuantumToolbox
using Random

# Global definition of the system
N = 10
a = kron(destroy(N), qeye(2))
σm = kron(qeye(N), sigmam())
σz = qeye(N) ⊗ sigmaz()

g = 0.01
ωc = 1
ωq = 0.99
γ = 0.1
nth = 0.001

# Jaynes-Cummings Hamiltonian
H = ωc * a' * a + ωq / 2 * σz + g * (a' * σm + a * σm')
ψ0 = kron(fock(N, 0), fock(2, 0))

e_ops = (a' * a, σz)
c_ops = (sqrt(γ * (1 + nth)) * a, sqrt(γ * nth) * a', sqrt(γ * (1 + nth)) * σm, sqrt(γ * nth) * σm')

sme_η = 0.7 # Efficiency of the homodyne detector for smesolve
c_ops_sme = ntuple(i -> sqrt(1 - sme_η) * c_ops[i], Val(length(c_ops)))
sc_ops_sme = ntuple(i -> sqrt(sme_η) * c_ops[i], Val(length(c_ops)))

# The following definition is to test the case of `sc_ops` as an `AbstractQuantumObject`
c_ops_sme2 = c_ops[2:end]
sc_ops_sme2 = c_ops[1]

ψ0_int = Qobj(round.(Int, real.(ψ0.data)), dims = ψ0.dims) # Used for testing the type inference

ψ_wrong = kron(fock(N - 1, 0), fock(2, 0))

rng = MersenneTwister(12)

# QobjEvo
ωd = 1.02
F = 0.05
coef1(p, t) = p.F * exp(1im * p.ωd * t)
coef2(p, t) = p.F * exp(-1im * p.ωd * t)
p = (F = F, ωd = ωd)
H_td = (H, (a, coef1), (a', coef2))
H_td2 = QobjEvo(H_td)
L_td = liouvillian(H_td2)

# time list and saveat
tlist = range(0, 10 / γ, 100)
saveat_idxs = 50:90
saveat = tlist[saveat_idxs]

# time list for testing exceptions
tlist1 = Float64[]
tlist2 = [0, 0.2, 0.1]
tlist3 = [0, 0.1, 0.1, 0.2]

# mesolve solution used for comparing results from mcsolve, ssesolve, and smesolve with mesolve
prob_me = mesolveProblem(H, ψ0, tlist, c_ops, e_ops = e_ops, progress_bar = Val(false))
sol_me = mesolve(prob_me)
end
