module QuantumToolboxCore

# standard Julia libraries
using SparseArrays
using LinearAlgebra
import LinearAlgebra: checksquare

# JuliaMath libraries
import TensorCore: TensorCore, tensor, ⊗

# other dependencies (in alphabetical order)
import Base: AbstractVecOrTuple
import FillArrays: Eye
import Random: AbstractRNG, default_rng
import StaticArraysCore: SVector, MVector

# SciML
import SciMLOperators:
    SciMLOperators,
    cache_operator,
    iscached,
    isconstant,
    AbstractSciMLOperator,
    MatrixOperator,
    ScalarOperator,
    ScaledOperator,
    AddedOperator,
    ComposedOperator,
    update_coefficients!,
    concretize

# Export functions from the other modules

## LinearAlgebra
export ishermitian, issymmetric, isposdef, dot, tr, svdvals, norm, normalize, normalize!, diag, Hermitian, Symmetric

## SparseArrays
export permute

## TensorCore
export tensor, ⊗

## SciMLOperators
export cache_operator, iscached, isconstant

## Basic utilities for QuantumToolbox libraries
include("settings.jl")
include("versioninfo.jl")
include("type_handle.jl")

## Linear Algebra
include("linalg/linalg.jl")
include("linalg/arnoldi.jl")
include("linalg/linear_maps.jl")

## Other physical related utilities
include("PhysicalConstants.jl")
include("physics_func.jl")

## Quantum Object
include("qobj/dimensions.jl")
include("qobj/energy_restricted.jl")
include("qobj/quantum_object_base.jl")
include("qobj/quantum_object.jl")
include("qobj/quantum_object_evo.jl")
include("qobj/boolean_functions.jl")
include("qobj/arithmetic_and_attributes.jl")
include("qobj/eigen.jl")
include("qobj/functions.jl")
include("qobj/states.jl")
include("qobj/operators.jl")
include("qobj/superoperators.jl")
include("qobj/synonyms.jl")

## Other functionalities for QuantumObject
include("entropy.jl")
include("metrics.jl")
include("negativity.jl")

## deprecated functions
include("deprecated.jl")

function __init__()
    # register QuantumToolbox library and its dependencies
    if (QuantumToolboxCore ∉ QT_LIBRARIES)
        # use pushfirst! so that main API libraries are at the front of the registry (for better display order in versioninfo)
        pushfirst!(QT_LIBRARIES, QuantumToolboxCore)

        # dependencies
        m_list = Module[SciMLOperators]
        foreach(m_list) do m
            (m ∉ DEP_PKGS) && push!(DEP_PKGS, m)
        end
    end
    return nothing
end

end
