module QuantumToolboxOperators

using LinearAlgebra
using SparseArrays

import SciMLOperators: AbstractSciMLOperator, AdjointOperator, IdentityOperator, ScaledOperator
import SciMLOperators: islinear, has_adjoint, concretize, cache_operator, iscached

export BosonicOperator, DestroyOperator, NumberOperator, DestroyPowerOperator, NormalOrderedOperator
export KroneckerOperator

include("bosonic.jl")
include("kronecker.jl")

end # module QuantumToolboxOperators
