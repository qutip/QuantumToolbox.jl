module QuantumToolboxCUDAExt

using QuantumToolbox
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR
import SparseArrays: SparseVector, SparseMatrixCSC

CuArray{T}(A::QuantumObject{Tq}) where {T,Tq<:Union{Vector,Matrix}} = QuantumObject(CuArray{T}(A.data), A.type, A.dims)
CuSparseVector{T}(A::QuantumObject{<:SparseVector}) where T = QuantumObject(CuSparseVector{T}(A.data), A.type, A.dims)
CuSparseMatrixCSC{T}(A::QuantumObject{<:SparseMatrixCSC}) where T = QuantumObject(CuSparseMatrixCSC{T}(A.data), A.type, A.dims)
CuSparseMatrixCSR{T}(A::QuantumObject{<:SparseMatrixCSC}) where T = QuantumObject(CuSparseMatrixCSR{T}(A.data), A.type, A.dims)

cu(A::QuantumObject; word_size::Int=32) = ((word_size == 64) || (word_size == 32)) ? cu(A, Val(word_size)) : throw(DomainError(word_size, "The word size should be 32 or 64."))
cu(A::QuantumObject{T}, word_size::TW) where {T<:Union{Vector,Matrix},TW<:Union{Val{32},Val{64}}} = CuArray{_change_eltype(eltype(A), word_size)}(A)
cu(A::QuantumObject{<:SparseVector}, word_size::TW) where TW<:Union{Val{32},Val{64}} = CuSparseVector{_change_eltype(eltype(A), word_size)}(A)
cu(A::QuantumObject{<:SparseMatrixCSC}, word_size::TW) where TW<:Union{Val{32},Val{64}} = CuSparseMatrixCSC{_change_eltype(eltype(A), word_size)}(A)

_change_eltype(::Type{T}, ::Val{64}) where {T<:Int} = Int64
_change_eltype(::Type{T}, ::Val{32}) where {T<:Int} = Int32
_change_eltype(::Type{T}, ::Val{64}) where {T<:AbstractFloat} = Float64
_change_eltype(::Type{T}, ::Val{32}) where {T<:AbstractFloat} = Float32
_change_eltype(::Type{Complex{T}}, ::Val{64}) where {T<:Union{Int,AbstractFloat}} = ComplexF64
_change_eltype(::Type{Complex{T}}, ::Val{32}) where {T<:Union{Int,AbstractFloat}} = ComplexF32

end