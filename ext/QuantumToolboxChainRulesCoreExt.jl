module QuantumToolboxChainRulesCoreExt

using LinearAlgebra
import QuantumToolboxCore: QuantumObject
import ChainRulesCore: rrule, NoTangent, Tangent

function rrule(::Type{QuantumObject}, data, type, dimensions)
    obj = QuantumObject(data, type, dimensions)
    f_pullback(Δobj) = (NoTangent(), Δobj.data, NoTangent(), NoTangent())
    f_pullback(Δobj_data::AbstractArray) = (NoTangent(), Δobj_data, NoTangent(), NoTangent())
    return obj, f_pullback
end

end
