module QuantumToolboxOperators

using LinearAlgebra
using SparseArrays

import SciMLOperators: AbstractSciMLOperator, AdjointOperator, IdentityOperator
import SciMLOperators: islinear, has_adjoint, concretize

export DestroyOperator, NumberOperator, DestroyPowerOperator

include("bosonic.jl")

end # module QuantumToolboxOperators
