#= 
All boolean functions for checking the data or type in `QuantumObject`
=#

export isket, isbra, isoper, isoperbra, isoperket, issuper

"""
    isbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`BraQuantumObject`](@ref) state.
"""
isbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: BraQuantumObject

"""
    isket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`KetQuantumObject`](@ref) state.
"""
isket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: KetQuantumObject

"""
    isoper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorQuantumObject`](@ref) state.
"""
isoper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorQuantumObject

"""
    isoperbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorBraQuantumObject`](@ref) state.
"""
isoperbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorBraQuantumObject

"""
    isoperket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorKetQuantumObject`](@ref) state.
"""
isoperket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorKetQuantumObject

"""
    issuper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`SuperOperatorQuantumObject`](@ref) state.
"""
issuper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: SuperOperatorQuantumObject

LinearAlgebra.ishermitian(A::QuantumObject{<:AbstractArray{T}}) where {T} = ishermitian(A.data)
LinearAlgebra.issymmetric(A::QuantumObject{<:AbstractArray{T}}) where {T} = issymmetric(A.data)
LinearAlgebra.isposdef(A::QuantumObject{<:AbstractArray{T}}) where {T} = isposdef(A.data)
