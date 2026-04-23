module QuantumToolboxUtils

using LinearAlgebra
using SparseArrays
import StaticArraysCore: SVector, MVector
import Random: AbstractRNG

import SciMLOperators: SciMLOperators, AbstractSciMLOperator

include("settings.jl")
include("versioninfo.jl")
include("type_handle.jl")
include("PhysicalConstants.jl")
include("misc.jl")

function __init__()
    _register_qt_library!(QuantumToolboxUtils)
    return nothing
end

end
