module QuantumToolboxUtils

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
    _register_qt_library!(QuantumToolboxUtils)
    return nothing
end

end
