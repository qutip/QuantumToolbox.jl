module QuantumToolbox

# Re-export:
#   1. StaticArraysCore.SVector for the type of dims
#   2. basic functions in LinearAlgebra and SparseArrays
#   3. some functions in SciMLOperators
import Reexport: @reexport
@reexport import StaticArraysCore: SVector
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport import SciMLOperators: cache_operator, iscached, isconstant

# other functions in LinearAlgebra
import LinearAlgebra: BlasReal, BlasInt, BlasFloat, BlasComplex, checksquare
import LinearAlgebra.BLAS: @blasfunc
import LinearAlgebra.LAPACK: hseqr!

# SciML packages (for QobjEvo, OrdinaryDiffEq, and LinearSolve)
import SciMLBase:
    solve,
    solve!,
    init,
    reinit!,
    remake,
    u_modified!,
    ODEFunction,
    ODEProblem,
    SDEProblem,
    EnsembleProblem,
    EnsembleSerial,
    EnsembleThreads,
    EnsembleSplitThreads,
    EnsembleDistributed,
    FullSpecialize,
    CallbackSet,
    ContinuousCallback,
    DiscreteCallback,
    AbstractSciMLProblem,
    AbstractODEIntegrator
import StochasticDiffEq: StochasticDiffEqAlgorithm, SRA1
import SciMLOperators:
    SciMLOperators,
    AbstractSciMLOperator,
    MatrixOperator,
    ScalarOperator,
    ScaledOperator,
    AddedOperator,
    IdentityOperator,
    update_coefficients!,
    concretize
import SciMLStructures: isscimlstructure, ismutablescimlstructure, hasportion, canonicalize, replace, replace!, Tunable
import LinearSolve: LinearProblem, SciMLLinearSolveAlgorithm, KrylovJL_MINRES, KrylovJL_GMRES
import DiffEqBase: get_tstops
import DiffEqCallbacks: PeriodicCallback, PresetTimeCallback, TerminateSteadyState
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm
import OrdinaryDiffEqTsit5: Tsit5
import DiffEqNoiseProcess: RealWienerProcess!

# other dependencies (in alphabetical order)
import ArrayInterface: allowed_getindex, allowed_setindex!
import Distributed: RemoteChannel
import FFTW: fft, fftshift
import Graphs: connected_components, DiGraph
import IncompleteLU: ilu
import Pkg
import Random: AbstractRNG, default_rng, seed!
import SpecialFunctions: loggamma
import StaticArraysCore: MVector

# Setting the number of threads to 1 allows
# to achieve better performances for more massive parallelizations
BLAS.set_num_threads(1)

# Utility
include("utilities.jl")
include("versioninfo.jl")
include("progress_bar.jl")
include("linear_maps.jl")

# Quantum Object
include("qobj/quantum_object_base.jl")
include("qobj/quantum_object.jl")
include("qobj/quantum_object_evo.jl")
include("qobj/boolean_functions.jl")
include("qobj/arithmetic_and_attributes.jl")
include("qobj/eigsolve.jl")
include("qobj/functions.jl")
include("qobj/states.jl")
include("qobj/operators.jl")
include("qobj/superoperators.jl")
include("qobj/synonyms.jl")

# time evolution
include("time_evolution/time_evo_parameters.jl")
include("time_evolution/time_evolution.jl")
include("time_evolution/callback_helpers.jl")
include("time_evolution/mesolve.jl")
include("time_evolution/lr_mesolve.jl")
include("time_evolution/sesolve.jl")
include("time_evolution/mcsolve.jl")
include("time_evolution/ssesolve.jl")
include("time_evolution/time_evolution_dynamical.jl")

# Others
include("permutation.jl")
include("correlations.jl")
include("wigner.jl")
include("spin_lattice.jl")
include("arnoldi.jl")
include("metrics.jl")
include("negativity.jl")
include("steadystate.jl")

# deprecated functions
include("deprecated.jl")

end
