module QuantumToolboxVisualization

## Re-export of QuantumToolboxCore
import Reexport: @reexport
@reexport using QuantumToolboxCore

import QuantumToolboxCore: makeVal

## dependencies (in alphabetical order)
import LaTeXStrings: @L_str
import LinearAlgebra: lmul!
import SparseArrays: AbstractSparseArray, findnz
import SpecialFunctions: loggamma

# Source files

## Some overloading with QuantumToolboxCore library
include("core_overload.jl")

## Visualization
include("bloch_sphere.jl")
include("fock_distribution.jl")
include("matrix.jl")
include("wigner.jl")

## deprecated functions
include("deprecated.jl")

function __init__()
    # register QuantumToolbox library and its dependencies
    if (QuantumToolboxVisualization ∉ QuantumToolboxCore.QT_LIBRARIES)
        # use pushfirst! so that main API libraries are at the front of the registry (for better display order in versioninfo)
        pushfirst!(QuantumToolboxCore.QT_LIBRARIES, QuantumToolboxVisualization)

        # so far, no need to add DEP_PKGS for this library
    end
    return nothing
end

end
