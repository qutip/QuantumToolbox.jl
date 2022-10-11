module QuPhys

using Reexport
using Distributed
@reexport using ProgressMeter
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using OrdinaryDiffEq
@reexport using DiffEqCallbacks

using OrdinaryDiffEq: OrdinaryDiffEqMutableCache

@inline function DiffEqBase.get_tmp_cache(integrator,
    alg::LinearExponential,
    cache::OrdinaryDiffEqMutableCache)
(cache.tmp,)
end

include("quantum_operators.jl")
include("general_functions.jl")
include("time_evolution.jl")

export spre, spost, sprepost, lindblad_dissipator
export fock, coherent
export destroy, create, sigmam, sigmap, eye, projection
export sinm, cosm
export expect
export wigner
export row_major_reshape, chop_op, gaussian, trunc_op, meshgrid
export ptrace, entropy_vn, entanglement
export mesolve, mcsolve, sesolve, liouvillian, liouvillian_floquet, steadystate, steadystate_floquet

end
