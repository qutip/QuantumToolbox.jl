module QuantumToolboxAMDGPUExt

using QuantumToolbox
using QuantumToolbox: makeVal, getVal
import QuantumToolbox: _sparse_similar, _convert_eltype_wordsize
import AMDGPU: roc, ROCArray, allowscalar
import AMDGPU.rocSPARSE: ROCSparseVector, ROCSparseMatrixCSC, ROCSparseMatrixCSR, AbstractROCSparseArray
import SparseArrays: SparseVector, SparseMatrixCSC, sparse
import AMDGPU.Adapt: adapt

allowscalar(false)

@doc raw"""
    ROCArray(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.ROCArray` for gpu calculations.
"""
ROCArray(A::QuantumObject) = QuantumObject(ROCArray(A.data), A.type, A.dimensions)

@doc raw"""
    ROCArray{T}(A::QuantumObject)

If `A.data` is a dense array, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.ROCArray` with element type `T` for gpu calculations.
"""
ROCArray{T}(A::QuantumObject) where {T} = QuantumObject(ROCArray{T}(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseVector(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseVector` for gpu calculations.
"""
ROCSparseVector(A::QuantumObject) = QuantumObject(ROCSparseVector(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseVector{T}(A::QuantumObject)

If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseVector` with element type `T` for gpu calculations.
"""
ROCSparseVector{T}(A::QuantumObject) where {T} = QuantumObject(ROCSparseVector{T}(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseMatrixCSC(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseMatrixCSC` for gpu calculations.
"""
ROCSparseMatrixCSC(A::QuantumObject) = QuantumObject(ROCSparseMatrixCSC(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseMatrixCSC{T}(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseMatrixCSC` with element type `T` for gpu calculations.
"""
ROCSparseMatrixCSC{T}(A::QuantumObject) where {T} = QuantumObject(ROCSparseMatrixCSC{T}(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseMatrixCSR` for gpu calculations.
"""
ROCSparseMatrixCSR(A::QuantumObject) = QuantumObject(ROCSparseMatrixCSR(A.data), A.type, A.dimensions)

@doc raw"""
    ROCSparseMatrixCSR(A::QuantumObject)

If `A.data` is in the type of `SparseMatrixCSC`, return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU.rocSPARSE.ROCSparseMatrixCSR` with element type `T` for gpu calculations.
"""
ROCSparseMatrixCSR{T}(A::QuantumObject) where {T} = QuantumObject(ROCSparseMatrixCSR{T}(A.data), A.type, A.dimensions)

@doc raw"""
    roc(A::QuantumObject; word_size::Int=64)

Return a new [`QuantumObject`](@ref) where `A.data` is in the type of `AMDGPU` arrays for gpu calculations.

# Arguments
- `A::QuantumObject`: The [`QuantumObject`](@ref)
- `word_size::Int`: The word size of the element type of `A`, can be either `32` or `64`. Default to `64`.
"""
function roc(A::QuantumObject; word_size::Union{Val,Int} = Val(64))
    _word_size = getVal(makeVal(word_size))

    ((_word_size == 64) || (_word_size == 32)) || throw(DomainError(_word_size, "The word size should be 32 or 64."))

    return roc(A, makeVal(word_size))
end
roc(A::QuantumObject, word_size::Union{Val{32},Val{64}}) =
    QuantumObject(adapt(ROCArray{_convert_eltype_wordsize(eltype(A), word_size)}, A.data), A.type, A.dimensions)
function roc(
    A::QuantumObject{ObjType,DimsType,<:SparseVector},
    word_size::Union{Val{32},Val{64}},
) where {ObjType<:QuantumObjectType,DimsType<:AbstractDimensions}
    return ROCSparseVector{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end
function roc(
    A::QuantumObject{ObjType,DimsType,<:SparseMatrixCSC},
    word_size::Union{Val{32},Val{64}},
) where {ObjType<:QuantumObjectType,DimsType<:AbstractDimensions}
    return ROCSparseMatrixCSC{_convert_eltype_wordsize(eltype(A), word_size)}(A)
end

QuantumToolbox.to_dense(A::MT) where {MT<:AbstractROCSparseArray} = ROCArray(A)

QuantumToolbox.to_dense(::Type{T1}, A::ROCArray{T2}) where {T1<:Number,T2<:Number} = ROCArray{T1}(A)
QuantumToolbox.to_dense(::Type{T}, A::AbstractROCSparseArray) where {T<:Number} = ROCArray{T}(A)

QuantumToolbox._sparse_similar(A::ROCSparseMatrixCSC, args...) = sparse(args..., fmt = :csc)
QuantumToolbox._sparse_similar(A::ROCSparseMatrixCSR, args...) = sparse(args..., fmt = :csr)

end
