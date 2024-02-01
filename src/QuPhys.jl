module QuPhys

using Reexport
using Distributed
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using OrdinaryDiffEq
@reexport using DiffEqCallbacks
using Graphs
using FFTW
using SpecialFunctions
using LinearSolve
using LinearMaps: LinearMap
using IncompleteLU

using LinearAlgebra: BlasFloat, BlasComplex

# Setting the number of threads to 1 allows
# to achieve better performances for more massive parallelizations
BLAS.set_num_threads(1)

include("quantum_object.jl")
include("quantum_operators.jl")
include("general_functions.jl")
include("time_evolution/time_evolution.jl")
include("time_evolution/mesolve.jl")
include("time_evolution/lr_mesolve.jl")
include("time_evolution/sesolve.jl")
include("time_evolution/mcsolve.jl")
include("time_evolution/time_evolution_dynamical.jl")
include("permutation.jl")
include("correlations.jl")
include("wigner.jl")
include("spin_lattice.jl")
include("arnoldi.jl")
include("eigsolve.jl")

export QuantumObject, Qobj, BraQuantumObject, KetQuantumObject, OperatorQuantumObject, SuperOperatorQuantumObject, TimeEvolutionSol
export isket, isbra, isoper, issuper, ket2dm
export spre, spost, sprepost, lindblad_dissipator
export fock, basis, coherent
export sigmam, sigmap, sigmax, sigmay, sigmaz
export destroy, create, eye, projection, rand_dm
export tensor, ⊗
export sinm, cosm
export expect
export WignerClenshaw, WignerLaguerre, wigner
export row_major_reshape, tidyup, tidyup!, gaussian, trunc_op, meshgrid, sparse_to_dense, dense_to_sparse
export get_data, mat2vec, vec2mat
export ptrace, entropy_vn, entanglement
export get_coherence, n_th
export dfd_mesolve, dsf_mesolve, dsf_mcsolve
export liouvillian, liouvillian_floquet, liouvillian_generalized, steadystate, steadystate_floquet
export LiouvillianDirectSolver, SteadyStateDirectSolver
export bdf, get_bdf_blocks
export FFTCorrelation, ExponentialSeries
export correlation_3op_2t, correlation_2op_2t, correlation_2op_1t, spectrum
export eigsolve, eigsolve_al
end
