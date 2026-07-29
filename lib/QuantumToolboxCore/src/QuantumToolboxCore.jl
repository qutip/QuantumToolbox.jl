module QuantumToolboxCore

using LinearAlgebra
using SparseArrays
import StaticArraysCore: SVector, MVector
import Random: AbstractRNG

import SciMLOperators: SciMLOperators, AbstractSciMLOperator

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
