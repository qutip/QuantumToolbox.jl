module QuantumToolbox

# Re-export:
#   1. basic functions in LinearAlgebra and SparseArrays
#   2. the solvers in ODE and LinearSolve
import Reexport: @reexport
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using OrdinaryDiffEq
@reexport using StochasticDiffEq
@reexport using LinearSolve

# other functions in LinearAlgebra
import LinearAlgebra: BlasReal, BlasInt, BlasFloat, BlasComplex, checksquare
import LinearAlgebra.BLAS: @blasfunc
if VERSION < v"1.10"
    import LinearAlgebra: chkstride1
    import LinearAlgebra.BLAS: libblastrampoline
    import LinearAlgebra.LAPACK: chklapackerror
    import Base: require_one_based_indexing
else
    import LinearAlgebra.LAPACK: hseqr!
end

# other dependencies (in alphabetical order)
import ArrayInterface: allowed_getindex, allowed_setindex!
import DiffEqCallbacks: DiscreteCallback, PeriodicCallback, PresetTimeCallback, TerminateSteadyState
import DiffEqNoiseProcess: RealWienerProcess, RealWienerProcess!
import FFTW: fft, fftshift
import Graphs: connected_components, DiGraph
import IncompleteLU: ilu
import LinearMaps: LinearMap
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm
import Pkg
import Random
import SpecialFunctions: loggamma

# Setting the number of threads to 1 allows
# to achieve better performances for more massive parallelizations
BLAS.set_num_threads(1)

# Utility
include("utilities.jl")
include("versioninfo.jl")
include("progress_bar.jl")

# Quantum Object
include("qobj/quantum_object.jl")
include("qobj/boolean_functions.jl")
include("qobj/arithmetic_and_attributes.jl")
include("qobj/eigsolve.jl")
include("qobj/functions.jl")
include("qobj/states.jl")
include("qobj/operators.jl")
include("qobj/superoperators.jl")
include("qobj/synonyms.jl")
include("qobj/operator_sum.jl")

# time evolution
include("time_evolution/time_evolution.jl")
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

end
