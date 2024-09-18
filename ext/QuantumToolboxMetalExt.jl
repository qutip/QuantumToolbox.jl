module QuantumToolboxMetalExt

using QuantumToolbox
import Metal: mtl, MtlArray

@doc raw"""
    MtlArray(A::QuantumObject)

If `A.data` is an arbitrary array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `Metal.MtlArray` for gpu calculations.

Note that this function will always change element type into `32`-bit (`Int32`, `Float32`, and `ComplexF32`).
"""
MtlArray(A::QuantumObject{<:AbstractArray{T}}) where {T<:Number} = QuantumObject(MtlArray(A.data), A.type, A.dims)
MtlArray(A::QuantumObject{<:AbstractArray{T}}) where {T<:Int64} = QuantumObject(MtlArray{Int32}(A.data), A.type, A.dims)
MtlArray(A::QuantumObject{<:AbstractArray{T}}) where {T<:Float64} =
    QuantumObject(MtlArray{Float32}(A.data), A.type, A.dims)
MtlArray(A::QuantumObject{<:AbstractArray{T}}) where {T<:ComplexF64} =
    QuantumObject(MtlArray{ComplexF32}(A.data), A.type, A.dims)

@doc raw"""
    MtlArray{T}(A::QuantumObject)

If `A.data` is an arbitrary array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `Metal.MtlArray` with element type `T` for gpu calculations.
"""
MtlArray{T}(A::QuantumObject{<:AbstractArray{Tq}}) where {T,Tq<:Number} =
    QuantumObject(MtlArray{T}(A.data), A.type, A.dims)

@doc raw"""
    mtl(A::QuantumObject)

Return a new [`QuantumObject`](@ref) where `A.data` is in the type of `Metal` arrays for gpu calculations.

Note that this function will always change element type into `32`-bit (`Int32`, `Float32`, and `ComplexF32`).
"""
mtl(A::QuantumObject{<:AbstractArray{T}}) where {T<:Int64} = QuantumObject(MtlArray{Int32}(A.data), A.type, A.dims)
mtl(A::QuantumObject{<:AbstractArray{T}}) where {T<:Float64} = QuantumObject(MtlArray{Float32}(A.data), A.type, A.dims)
mtl(A::QuantumObject{<:AbstractArray{T}}) where {T<:ComplexF64} =
    QuantumObject(MtlArray{ComplexF32}(A.data), A.type, A.dims)

## TODO: Remove the following part if Metal.jl support `sparse`
import QuantumToolbox: _spre, _spost, _sprepost
_spre(A::MtlArray, Id::AbstractMatrix) = kron(Id, A)
_spost(B::MtlArray, Id::AbstractMatrix) = kron(transpose(B), Id)
_sprepost(A::MtlArray, B::MtlArray) = kron(transpose(B), A)

## TODO: Remove the following part if Metal.jl support `kron`
import LinearAlgebra: kron
LinearAlgebra.kron(A::Diagonal, B::MtlArray{T}) where {T} = MtlArray{T}(LinearAlgebra.kron(mtl(A), B))
LinearAlgebra.kron(A::MtlArray{T}, B::Diagonal) where {T} = MtlArray{T}(LinearAlgebra.kron(A, mtl(B)))
LinearAlgebra.kron(A::Transpose{T1,<:MtlArray}, B::MtlArray{T2}) where {T1,T2} = MtlArray(LinearAlgebra.kron(A, B))
LinearAlgebra.kron(A::MtlArray{T1}, B::Transpose{T2,<:MtlArray}) where {T1,T2} = MtlArray(LinearAlgebra.kron(A, B))
LinearAlgebra.kron(A::Transpose{T1,<:MtlArray}, B::Transpose{T2,<:MtlArray}) where {T1,T2} =
    MtlArray(LinearAlgebra.kron(A, B))
LinearAlgebra.kron(A::MtlArray{T1}, B::MtlArray{T2}) where {T1,T2} = MtlArray(LinearAlgebra.kron(A, B))

end
