module QuPhys

using LinearAlgebra
using SparseArrays

include("quantum_operators.jl")
include("general_functions.jl")

export spre, spost, sprepost, lindblad_dissipator
export destroy, eye, fock, projection
export sinm, cosm
export expect
export chop_op, gaussian, gaussian_derivative, trunc_op

end
