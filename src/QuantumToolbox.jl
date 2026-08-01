module QuantumToolbox

## Standard Julia libraries
using LinearAlgebra
using SparseArrays

import Distributed: RemoteChannel
import Random: AbstractRNG, default_rng
import Statistics: mean, std

## Re-export of QuantumToolbox libraries
import Reexport: @reexport
@reexport using QuantumToolboxCore

## internal functions of QuantumToolbox libraries
import QuantumToolboxCore:
    position, # since we don't export it (conflict with Base.position)
    momentum, # since we don't export it (just to align with position)
    getVal,
    makeVal,
    get_typename_wrapper,
    isendomorphic,
    promote_op_type,
    check_dimensions,
    check_mul_dimensions,
    dimensions_to_dims,
    to_sparse_if_needed,
    _float_type,
    _complex_float_type,
    _non_endomorphic_dims_error,
    _dense_similar,
    _sparse_similar,
    _spre,
    _spost,
    _sprepost,
    _ptrace_oper

## SciML packages (for QobjEvo, OrdinaryDiffEq, and LinearSolve)
import SciMLBase:
    SciMLBase,
    solve,
    solve!,
    init,
    reinit!,
    remake,
    derivative_discontinuity!,
    NullParameters,
    LinearProblem,
    ODEFunction,
    SDEFunction,
    ODEProblem,
    SDEProblem,
    EnsembleProblem,
    EnsembleContext,
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
import StochasticDiffEqHighOrder: SRA2, SRIW1
import SciMLOperators:
    cache_operator,
    isconstant,
    SciMLOperators,
    AbstractSciMLOperator,
    MatrixOperator,
    ScalarOperator,
    AddedOperator,
    IdentityOperator
import LinearSolve:
    LinearSolve, SciMLLinearSolveAlgorithm, KrylovJL_MINRES, KrylovJL_GMRES, UMFPACKFactorization, LUFactorization, OperatorAssumptions
import DiffEqCallbacks: PeriodicCallback, FunctionCallingCallback, FunctionCallingAffect, TerminateSteadyState
import OrdinaryDiffEqCore
import OrdinaryDiffEqVerner: Vern7
import OrdinaryDiffEqLowOrderRK: DP5
import DiffEqNoiseProcess: RealWienerProcess!, RealWienerProcess

## other dependencies (in alphabetical order)
import ArrayInterface: allowed_setindex!
import FFTW: fft, ifft, fftfreq, fftshift
import FillArrays: Eye
import Graphs: connected_components, DiGraph
import IncompleteLU: ilu
import LaTeXStrings: @L_str
import ProgressMeter: Progress, next!
import SpecialFunctions: loggamma
import StaticArraysCore: SVector, MVector

# Source files

## Some overloading with QuantumToolboxCore library
include("core_overload.jl")

## some functions for Quantum Object
include("eigsolve.jl")
include("block_diagonal_form.jl")

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
include("time_evolution/liouvillian_dressed_nonsecular.jl")
include("time_evolution/time_evolution_dynamical.jl")

## Other functionalities
include("correlations.jl")
include("wigner.jl")
include("spin_lattice.jl")
include("steadystate.jl")
include("spectrum.jl")

## Visualization
include("visualization/bloch_sphere.jl")
include("visualization/fock_distribution.jl")
include("visualization/matrix.jl")
include("visualization/wigner.jl")

## deprecated functions
include("deprecated.jl")

function __init__()
    # register QuantumToolbox library and its dependencies
    if (QuantumToolbox ∉ QuantumToolboxCore.QT_LIBRARIES)
        # use pushfirst! so that main API libraries are at the front of the registry (for better display order in versioninfo)
        pushfirst!(QuantumToolboxCore.QT_LIBRARIES, QuantumToolbox)

        # dependencies
        m_list = Module[SciMLBase, SciMLOperators, OrdinaryDiffEqCore, LinearSolve]
        foreach(m_list) do m
            (m ∉ QuantumToolboxCore.DEP_PKGS) && push!(QuantumToolboxCore.DEP_PKGS, m)
        end
    end

    return nothing
end

end
