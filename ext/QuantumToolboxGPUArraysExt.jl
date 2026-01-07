module QuantumToolboxGPUArraysExt

using QuantumToolbox

import GPUArrays: AbstractGPUArray, AbstractCuSparseArray, dense_array_type
import KernelAbstractions
import KernelAbstractions: @kernel, @Const, @index, get_backend, synchronize

QuantumToolbox.to_dense(::Type{T1}, A::AbstractGPUArray{T2}) where {T1<:Number,T2<:Number} = dense_array_type(A){T1}(A)
QuantumToolbox.to_dense(A::AbstractCuSparseArray) = dense_array_type(A)(A)
QuantumToolbox.to_dense(::Type{T}, A::AbstractCuSparseArray) where {T<:Number} = dense_array_type(A){T}(A)

@kernel function tr_kernel!(B, @Const(A))
    # i, j, k = @index(Global, NTuple)
    # Atomix.@atomic B[i, j] += A[i, j, k, k] # TODO: use Atomix when it will support Complex types

    i, j = @index(Global, NTuple)
    @inbounds B[i, j] = 0
    @inbounds for k in 1:size(A, 3)
        B[i, j] += A[i, j, k, k]
    end
end

function QuantumToolbox._map_trace(A::AbstractGPUArray{T, 4}) where {T}
    B = similar(A, size(A, 1), size(A, 2))
    fill!(B, 0)

    backend = get_backend(A)
    kernel! = tr_kernel!(backend)

    kernel!(B, A, ndrange = size(A)[1:2])
    KernelAbstractions.synchronize(backend)

    return B
end

end
