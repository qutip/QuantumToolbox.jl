module QuantumToolboxChainRulesCore

using LinearAlgebra
import QuantumToolbox: QuantumObject
using ChainRulesCore

function ChainRulesCore.rrule(::Type{QuantumObject}, data, type, dimensions)
    obj = QuantumObject(data, type, dimensions)
    f_pullback(Δobj) = (NoTangent(), Δobj.data, NoTangent(), NoTangent())
    f_pullback(Δobj_data::AbstractArray) = (NoTangent(), Δobj_data, NoTangent(), NoTangent())
    return obj, f_pullback
end

end
