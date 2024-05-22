module QuantumToolboxCUDAExt

using QuantumToolbox
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR
import SparseArrays: SparseVector, SparseMatrixCSC

@doc raw"""
    CuArray(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CuArray` for gpu calculations.
"""
CuArray(A::QuantumObject{Tq}) where {Tq<:Union{Vector,Matrix}} = QuantumObject(CuArray(A.data), A.type, A.dims)

@doc raw"""
    CuArray{T}(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CuArray` with element type `T` for gpu calculations.
"""
CuArray{T}(A::QuantumObject{Tq}) where {T,Tq<:Union{Vector,Matrix}} = QuantumObject(CuArray{T}(A.data), A.type, A.dims)

@doc raw"""
    CuSparseVector(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseVector` for gpu calculations.
"""
CuSparseVector(A::QuantumObject{<:SparseVector}) = QuantumObject(CuSparseVector(A.data), A.type, A.dims)

@doc raw"""
    CuSparseVector{T}(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseVector` with element type `T` for gpu calculations.
"""
CuSparseVector{T}(A::QuantumObject{<:SparseVector}) where {T} = QuantumObject(CuSparseVector{T}(A.data), A.type, A.dims)

@doc raw"""
    CuSparseMatrixCSC(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSC` for gpu calculations.
"""
CuSparseMatrixCSC(A::QuantumObject{<:SparseMatrixCSC}) = QuantumObject(CuSparseMatrixCSC(A.data), A.type, A.dims)

@doc raw"""
    CuSparseMatrixCSC{T}(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSC` with element type `T` for gpu calculations.
"""
CuSparseMatrixCSC{T}(A::QuantumObject{<:SparseMatrixCSC}) where {T} =
    QuantumObject(CuSparseMatrixCSC{T}(A.data), A.type, A.dims)

@doc raw"""
    CuSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` for gpu calculations.
"""
CuSparseMatrixCSR(A::QuantumObject{<:SparseMatrixCSC}) = QuantumObject(CuSparseMatrixCSR(A.data), A.type, A.dims)

@doc raw"""
    CuSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` with element type `T` for gpu calculations.
"""
CuSparseMatrixCSR{T}(A::QuantumObject{<:SparseMatrixCSC}) where {T} =
    QuantumObject(CuSparseMatrixCSR{T}(A.data), A.type, A.dims)

@doc raw"""
    cu(A::QuantumObject; word_size::Int=32)

Return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA` arrays for gpu calculations.

# Arguments
- `A::QuantumObject`: The [`QuantumObject`](@ref)
- `word_size::Int`: The word size of the element type of `A`, can be either `32` or `64`. Default to `32`.
"""
cu(A::QuantumObject; word_size::Int = 32) =
    ((word_size == 64) || (word_size == 32)) ? cu(A, Val(word_size)) :
    throw(DomainError(word_size, "The word size should be 32 or 64."))
cu(A::QuantumObject{T}, word_size::TW) where {T<:Union{Vector,Matrix},TW<:Union{Val{32},Val{64}}} =
    CuArray{_change_eltype(eltype(A), word_size)}(A)
cu(A::QuantumObject{<:SparseVector}, word_size::TW) where {TW<:Union{Val{32},Val{64}}} =
    CuSparseVector{_change_eltype(eltype(A), word_size)}(A)
cu(A::QuantumObject{<:SparseMatrixCSC}, word_size::TW) where {TW<:Union{Val{32},Val{64}}} =
    CuSparseMatrixCSC{_change_eltype(eltype(A), word_size)}(A)

_change_eltype(::Type{T}, ::Val{64}) where {T<:Int} = Int64
_change_eltype(::Type{T}, ::Val{32}) where {T<:Int} = Int32
_change_eltype(::Type{T}, ::Val{64}) where {T<:AbstractFloat} = Float64
_change_eltype(::Type{T}, ::Val{32}) where {T<:AbstractFloat} = Float32
_change_eltype(::Type{Complex{T}}, ::Val{64}) where {T<:Union{Int,AbstractFloat}} = ComplexF64
_change_eltype(::Type{Complex{T}}, ::Val{32}) where {T<:Union{Int,AbstractFloat}} = ComplexF32

end
