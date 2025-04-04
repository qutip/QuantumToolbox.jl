module QuantumToolboxCUDAExt

using QuantumToolbox
using QuantumToolbox: makeVal, getVal
import QuantumToolbox: _sparse_similar, _convert_eltype_wordsize
import CUDA: cu, CuArray, allowscalar
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseArray
import SparseArrays: SparseVector, SparseMatrixCSC
import CUDA.Adapt: adapt

allowscalar(false)

@doc raw"""
    CuArray(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CuArray` for gpu calculations.
"""
CuArray(A::QuantumObject) = QuantumObject(CuArray(A.data), A.type, A.dimensions)

@doc raw"""
    CuArray{T}(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CuArray` with element type `T` for gpu calculations.
"""
CuArray{T}(A::QuantumObject) where {T} = QuantumObject(CuArray{T}(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseVector(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseVector` for gpu calculations.
"""
CuSparseVector(A::QuantumObject) = QuantumObject(CuSparseVector(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseVector{T}(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseVector` with element type `T` for gpu calculations.
"""
CuSparseVector{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseVector{T}(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseMatrixCSC(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSC` for gpu calculations.
"""
CuSparseMatrixCSC(A::QuantumObject) = QuantumObject(CuSparseMatrixCSC(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseMatrixCSC{T}(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSC` with element type `T` for gpu calculations.
"""
CuSparseMatrixCSC{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseMatrixCSC{T}(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` for gpu calculations.
"""
CuSparseMatrixCSR(A::QuantumObject) = QuantumObject(CuSparseMatrixCSR(A.data), A.type, A.dimensions)

@doc raw"""
    CuSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` with element type `T` for gpu calculations.
"""
CuSparseMatrixCSR{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseMatrixCSR{T}(A.data), A.type, A.dimensions)

@doc raw"""
    cu(A::QuantumObject; word_size::Int=64)

Return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA` arrays for gpu calculations.

# Arguments
- `A::QuantumObject`: The [`QuantumObject`](@ref)
- `word_size::Int`: The word size of the element type of `A`, can be either `32` or `64`. Default to `64`.
"""
function cu(A::QuantumObject; word_size::Union{Val,Int} = Val(64))
    _word_size = getVal(makeVal(word_size))

    ((_word_size == 64) || (_word_size == 32)) || throw(DomainError(_word_size, "The word size should be 32 or 64."))

    return cu(A, makeVal(word_size))
end
cu(A::QuantumObject, word_size::Union{Val{32},Val{64}}) = QuantumObject(adapt(CuArray{_convert_eltype_wordsize(eltype(A), word_size)}, A.data), A.type, A.dimensions)
function cu(
    A::QuantumObject{ObjType,DimsType,<:SparseVector},
    word_size::Union{Val{32},Val{64}},
) where {ObjType<:QuantumObjectType,DimsType<:AbstractDimensions}
    return CuSparseVector{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end
function cu(
    A::QuantumObject{ObjType,DimsType,<:SparseMatrixCSC},
    word_size::Union{Val{32},Val{64}},
) where {ObjType<:QuantumObjectType,DimsType<:AbstractDimensions}
    return CuSparseMatrixCSC{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end

QuantumToolbox.to_dense(A::MT) where {MT<:AbstractCuSparseArray} = CuArray(A)

QuantumToolbox.to_dense(::Type{T1}, A::CuArray{T2}) where {T1<:Number,T2<:Number} = CuArray{T1}(A)
QuantumToolbox.to_dense(::Type{T}, A::AbstractCuSparseArray) where {T<:Number} = CuArray{T}(A)

QuantumToolbox._sparse_similar(A::CuSparseMatrixCSC, args...) = sparse(args..., fmt = :csc)
QuantumToolbox._sparse_similar(A::CuSparseMatrixCSR, args...) = sparse(args..., fmt = :csr)

end
