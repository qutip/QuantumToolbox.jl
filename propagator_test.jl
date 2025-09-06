using QuantumToolbox
using LinearAlgebra

using Revise

# p_se = propagator(QobjEvo(sigmax()))
# p_me = propagator(sigmax(), [sigmaz()])

# println(p_se(1.0))
# println("-------------------")
# println(p_me(2, 0.1))

sx = sigmax()
sy = sigmay()
sz = sigmaz()
I = one(sz)

H = sz + QobjEvo((0.1*sx, (p, t) -> sin(2*t)))

p_se = propagator(H, progress_bar = false)

c_ops = [0.1*sz]
p_me = propagator(H, c_ops; progress_bar = false)

sesolve_res = []
mesolve_res = []

t0 = 0.0
t = 10.0
for i in 0:1
    psi = fock(2, i)
    push!(sesolve_res, sesolve(H, psi, [t0, t], progress_bar = false))
end

sup_psi_base = mat2vec(fock(2, 0)*fock(2, 0)')
for i in 1:4
    sup_psi = 0 * sup_psi_base
    sup_psi[i] = 1
    push!(mesolve_res, mesolve(H, vec2mat(sup_psi), [t0, t], c_ops, progress_bar = false))
end

U_se = p_se(t)
U_me = p_me(t)

U_sesolve = QuantumObject(hcat([sesolve_res[i].states[end].data for i in 1:2]...))
U_mesolve = QuantumObject(hcat([mat2vec(mesolve_res[i].states[end]).data for i in 1:4]...), type = SuperOperator())

println(norm(U_se-U_sesolve)<1e-8)
println(norm(U_me-U_mesolve)<1e-8)
