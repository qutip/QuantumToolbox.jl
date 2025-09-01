Pkg.activate(".")
include("src/QuantumToolbox.jl")
using .QuantumToolbox
using LinearAlgebra

using Revise


p_se = propagator(sigmax())
p_me = propagator(sigmax(), [sigmaz()])

println(p_se(1.0))
println("-------------------")
println(p_me(2, t0 = 1))