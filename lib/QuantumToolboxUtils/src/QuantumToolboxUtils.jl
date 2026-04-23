module QuantumToolboxUtils

using LinearAlgebra
using SparseArrays
import StaticArraysCore: SVector, MVector
import Random: AbstractRNG

import SciMLOperators: SciMLOperators, AbstractSciMLOperator

include("settings.jl")
include("versioninfo.jl")

include("type_handle.jl")
include("linalg.jl")
include("linear_maps.jl")
include("PhysicalConstants.jl")
include("physics_func.jl")

## deprecated functions
include("deprecated.jl")

function __init__()
    _register_qt_library!(QuantumToolboxUtils)
    return nothing
end

end
