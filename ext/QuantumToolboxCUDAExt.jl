module QuantumToolboxCUDAExt

# QuantumToolbox
using QuantumToolbox
using QuantumToolbox: makeVal, getVal
import QuantumToolbox: _sparse_similar, _convert_eltype_wordsize

# CUDA libraries
import CUDACore: CUDACore, cu, CuArray, allowscalar
import CUDACore.Adapt: adapt
import cuSPARSE: cuSPARSE, CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR

# other imports
import SparseArrays: SparseVector, SparseMatrixCSC, sparse

CUDACore.allowscalar(false)

@doc raw"""
    CUDACore.CuArray(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDACore.CuArray` for gpu calculations.
"""
CUDACore.CuArray(A::QuantumObject) = QuantumObject(CuArray(A.data), A.type, A.dimensions)

@doc raw"""
    CUDACore.CuArray{T}(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDACore.CuArray` with element type `T` for gpu calculations.
"""
CUDACore.CuArray{T}(A::QuantumObject) where {T} = QuantumObject(CuArray{T}(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseVector(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseVector` for gpu calculations.
"""
cuSPARSE.CuSparseVector(A::QuantumObject) = QuantumObject(CuSparseVector(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseVector{T}(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseVector` with element type `T` for gpu calculations.
"""
cuSPARSE.CuSparseVector{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseVector{T}(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseMatrixCSC(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseMatrixCSC` for gpu calculations.
"""
cuSPARSE.CuSparseMatrixCSC(A::QuantumObject) = QuantumObject(CuSparseMatrixCSC(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseMatrixCSC{T}(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseMatrixCSC` with element type `T` for gpu calculations.
"""
cuSPARSE.CuSparseMatrixCSC{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseMatrixCSC{T}(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseMatrixCSR` for gpu calculations.
"""
cuSPARSE.CuSparseMatrixCSR(A::QuantumObject) = QuantumObject(CuSparseMatrixCSR(A.data), A.type, A.dimensions)

@doc raw"""
    cuSPARSE.CuSparseMatrixCSR{T}(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `cuSPARSE.CuSparseMatrixCSR` with element type `T` for gpu calculations.
"""
cuSPARSE.CuSparseMatrixCSR{T}(A::QuantumObject) where {T} = QuantumObject(CuSparseMatrixCSR{T}(A.data), A.type, A.dimensions)

@doc raw"""
    CUDACore.cu(A::QuantumObject; word_size::Int=64)

Return a new [`QuantumObject`](@ref) where `A.data` is in the type of `CUDA` (dense or sparse) arrays for gpu calculations.

# Arguments
- `A::QuantumObject`: The [`QuantumObject`](@ref)
- `word_size::Int`: The word size of the element type of `A`, can be either `32` or `64`. Default to `64`.
"""
function CUDACore.cu(A::QuantumObject; word_size::Union{Val, Int} = Val(64))
    _word_size = getVal(makeVal(word_size))

    ((_word_size == 64) || (_word_size == 32)) || throw(DomainError(_word_size, "The word size should be 32 or 64."))

    return cu(A, makeVal(word_size))
end
CUDACore.cu(A::QuantumObject, word_size::Union{Val{32}, Val{64}}) =
    QuantumObject(adapt(CuArray{_convert_eltype_wordsize(eltype(A), word_size)}, A.data), A.type, A.dimensions)
function CUDACore.cu(
        A::QuantumObject{ObjType, DimsType, <:SparseVector},
        word_size::Union{Val{32}, Val{64}},
    ) where {ObjType <: QuantumObjectType, DimsType <: Dimensions}
    return CuSparseVector{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end
function CUDACore.cu(
        A::QuantumObject{ObjType, DimsType, <:SparseMatrixCSC},
        word_size::Union{Val{32}, Val{64}},
    ) where {ObjType <: QuantumObjectType, DimsType <: Dimensions}
    return CuSparseMatrixCSC{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end

QuantumToolbox._sparse_similar(A::CuSparseMatrixCSC, args...) = sparse(args..., fmt = :csc)
QuantumToolbox._sparse_similar(A::CuSparseMatrixCSR, args...) = sparse(args..., fmt = :csr)

end
