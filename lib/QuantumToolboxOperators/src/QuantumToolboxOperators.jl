module QuantumToolboxOperators

using LinearAlgebra
using SparseArrays

import SciMLOperators
import SciMLOperators: AbstractSciMLOperator, AdjointOperator, IdentityOperator, ScaledOperator
import SciMLOperators: islinear, has_adjoint, concretize, cache_operator, iscached

export BosonicOperator, DestroyOperator, NumberOperator, DestroyPowerOperator, NormalOrderedOperator
export LocalTensorProductOperator

include("bosonic.jl")
include("tensor_product.jl")

end
