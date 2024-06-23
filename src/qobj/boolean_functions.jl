#= 
All boolean functions for checking the data or type in `QuantumObject`
=#

export isket, isbra, isoper, isoperbra, isoperket, issuper

@doc raw"""
    isbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`BraQuantumObject`](@ref).
"""
isbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: BraQuantumObject

@doc raw"""
    isket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`KetQuantumObject`](@ref).
"""
isket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: KetQuantumObject

@doc raw"""
    isoper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorQuantumObject`](@ref).
"""
isoper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorQuantumObject

@doc raw"""
    isoperbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorBraQuantumObject`](@ref).
"""
isoperbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorBraQuantumObject

@doc raw"""
    isoperket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorKetQuantumObject`](@ref).
"""
isoperket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorKetQuantumObject

@doc raw"""
    issuper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`SuperOperatorQuantumObject`](@ref).
"""
issuper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: SuperOperatorQuantumObject

@doc raw"""
    ishermitian(A::QuantumObject)

Test whether the [`QuantumObject`](@ref) is Hermitian.
"""
LinearAlgebra.ishermitian(A::QuantumObject{<:AbstractArray{T}}) where {T} = ishermitian(A.data)

@doc raw"""
    issymmetric(A::QuantumObject)

Test whether the [`QuantumObject`](@ref) is symmetric.
"""
LinearAlgebra.issymmetric(A::QuantumObject{<:AbstractArray{T}}) where {T} = issymmetric(A.data)

@doc raw"""
    isposdef(A::QuantumObject)

Test whether the [`QuantumObject`](@ref) is positive definite (and Hermitian) by trying to perform a Cholesky factorization of `A`.
"""
LinearAlgebra.isposdef(A::QuantumObject{<:AbstractArray{T}}) where {T} = isposdef(A.data)
