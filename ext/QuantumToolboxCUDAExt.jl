module QuantumToolboxCUDAExt

using QuantumToolbox
import Sys
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR
import SparseArrays: SparseVector, SparseMatrixCSC

cu(A::QuantumObject; word_size::Int=Sys.WORD_SIZE) = ((word_size == 64) || (word_size == 32)) ? cu(A, Val(word_size)) : throw(DomainError(word_size, "The word size should be 32 or 64."))

cu(A::QuantumObject{T}, word_size::TW) where 
{T<:Union{Vector,Matrix},TW<:Union{Val{32},Val{64}}} = QuantumObject(CuArray{_change_eltype(eltype(A), word_size)}(A.data), A.type, A.dims)

function cu(A::QuantumObject{<:SparseVector}, word_size::TW) where TW<:Union{Val{32},Val{64}}
    A_cpu  = A.data
    Ti     = _change_eltype(Int, word_size)
    Tv     = _change_eltype(eltype(A_cpu), word_size)
    rowval = CuArray{Ti}(A_cpu.nzind)
    nzval  = CuArray{Tv}(A_cpu.nzval)
    return QuantumObject(CuSparseVector{Tv, Ti}(rowval, nzval, A_cpu.n), A.type, A.dims)
end

function cu(A::QuantumObject{<:SparseMatrixCSC}, word_size::TW) where TW<:Union{Val{32},Val{64}}
    A_cpu  = A.data
    Ti     = _change_eltype(Int, word_size)
    Tv     = _change_eltype(eltype(A_cpu), word_size)
    colptr = CuArray{Ti}(A_cpu.colptr)
    rowval = CuArray{Ti}(A_cpu.rowval)
    nzval  = CuArray{Tv}(A_cpu.nzval)
    return QuantumObject(CuSparseMatrixCSC{Tv, Ti}(colptr, rowval, nzval, size(A_cpu)), A.type, A.dims)
end

_change_eltype(::Type{T}, ::Val{64}) where {T<:Int} = Int64
_change_eltype(::Type{T}, ::Val{32}) where {T<:Int} = Int32
_change_eltype(::Type{T}, ::Val{64}) where {T<:AbstractFloat} = Float64
_change_eltype(::Type{T}, ::Val{32}) where {T<:AbstractFloat} = Float32
_change_eltype(::Type{Complex{T}}, ::Val{64}) where {T<:Union{Int,AbstractFloat}} = ComplexF64
_change_eltype(::Type{Complex{T}}, ::Val{32}) where {T<:Union{Int,AbstractFloat}} = ComplexF32

end