module QuantumToolbox

## Standard Julia libraries
using LinearAlgebra
using SparseArrays

import Distributed: RemoteChannel
import LinearAlgebra: checksquare
import Pkg
import Random: AbstractRNG, default_rng, seed!
import Statistics: mean, std

## SciML packages (for QobjEvo, OrdinaryDiffEq, and LinearSolve)
import SciMLBase:
    solve,
    solve!,
    init,
    reinit!,
    remake,
    u_modified!,
    NullParameters,
    LinearProblem,
    ODEFunction,
    SDEFunction,
    ODEProblem,
    SDEProblem,
    EnsembleProblem,
    EnsembleAlgorithm,
    EnsembleSerial,
    EnsembleThreads,
    EnsembleSplitThreads,
    EnsembleDistributed,
    FullSpecialize,
    CallbackSet,
    ContinuousCallback,
    DiscreteCallback,
    AbstractSciMLProblem,
    AbstractODEIntegrator,
    AbstractODEAlgorithm,
    AbstractODESolution,
    AbstractSDEAlgorithm
import StochasticDiffEq: SRA2, SRIW1
import SciMLOperators:
    cache_operator,
    iscached,
    isconstant,
    SciMLOperators,
    AbstractSciMLOperator,
    MatrixOperator,
    ScalarOperator,
    ScaledOperator,
    AddedOperator,
    ComposedOperator,
    IdentityOperator,
    update_coefficients!,
    concretize
import LinearSolve:
    SciMLLinearSolveAlgorithm, KrylovJL_MINRES, KrylovJL_GMRES, UMFPACKFactorization, LUFactorization, OperatorAssumptions
import DiffEqCallbacks: PeriodicCallback, FunctionCallingCallback, FunctionCallingAffect, TerminateSteadyState
import OrdinaryDiffEqVerner: Vern7
import OrdinaryDiffEqLowOrderRK: DP5
import DiffEqNoiseProcess: RealWienerProcess!, RealWienerProcess

## other dependencies (in alphabetical order)
import ArrayInterface: allowed_getindex, allowed_setindex!
import FFTW: fft, ifft, fftfreq, fftshift
import FillArrays: Eye
import Graphs: connected_components, DiGraph
import IncompleteLU: ilu
import LaTeXStrings: @L_str
import ProgressMeter: Progress, next!
import SpecialFunctions: loggamma
import StaticArraysCore: SVector, MVector

# Export functions from the other modules

## LinearAlgebra
export ishermitian, issymmetric, isposdef, dot, tr, svdvals, norm, normalize, normalize!, diag, Hermitian, Symmetric

## SparseArrays
export permute

## SciMLOperators
export cache_operator, iscached, isconstant

# Source files

## Utility
include("settings.jl")
include("utilities.jl")
include("versioninfo.jl")
include("linear_maps.jl")

## Quantum Object
include("qobj/dimensions.jl")
include("qobj/energy_restricted.jl")
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
include("qobj/block_diagonal_form.jl")

## time evolution
include("time_evolution/time_evolution.jl")
include("time_evolution/callback_helpers/callback_helpers.jl")
include("time_evolution/callback_helpers/sesolve_callback_helpers.jl")
include("time_evolution/callback_helpers/mesolve_callback_helpers.jl")
include("time_evolution/callback_helpers/mcsolve_callback_helpers.jl")
include("time_evolution/callback_helpers/ssesolve_callback_helpers.jl")
include("time_evolution/callback_helpers/smesolve_callback_helpers.jl")
include("time_evolution/mesolve.jl")
include("time_evolution/brmesolve.jl")
include("time_evolution/lr_mesolve.jl")
include("time_evolution/sesolve.jl")
include("time_evolution/mcsolve.jl")
include("time_evolution/ssesolve.jl")
include("time_evolution/smesolve.jl")
include("time_evolution/time_evolution_dynamical.jl")

## Other functionalities
include("correlations.jl")
include("wigner.jl")
include("spin_lattice.jl")
include("arnoldi.jl")
include("entropy.jl")
include("metrics.jl")
include("negativity.jl")
include("steadystate.jl")
include("spectrum.jl")

## Visualization
include("visualization/bloch_sphere.jl")
include("visualization/fock_distribution.jl")
include("visualization/matrix.jl")
include("visualization/wigner.jl")

## deprecated functions
include("deprecated.jl")

end
