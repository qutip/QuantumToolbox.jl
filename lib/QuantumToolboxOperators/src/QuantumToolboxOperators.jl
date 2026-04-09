module QuantumToolboxOperators

using LinearAlgebra
using SparseArrays

import SciMLOperators: AbstractSciMLOperator, AdjointOperator, IdentityOperator, ScaledOperator
import SciMLOperators: islinear, has_adjoint, concretize

export BosonicOperator, DestroyOperator, NumberOperator, DestroyPowerOperator, NormalOrderedOperator

include("bosonic.jl")

end # module QuantumToolboxOperators
