module QuantumToolboxGPUArraysExt

using QuantumToolbox

import GPUArrays: AbstractGPUArray
import KernelAbstractions
import KernelAbstractions: @kernel, @Const, @index, get_backend, synchronize

@kernel function tr_kernel!(B, @Const(A))
    # i, j, k = @index(Global, NTuple)
    # Atomix.@atomic B[i, j] += A[i, j, k, k] # TODO: use Atomix when it will support Complex types

    i, j = @index(Global, NTuple)
    @inbounds B[i, j] = 0
    @inbounds for k in 1:size(A, 3)
        B[i, j] += A[i, j, k, k]
    end
end

function QuantumToolbox._map_trace(A::AbstractGPUArray{T,4}) where {T}
    B = similar(A, size(A, 1), size(A, 2))
    fill!(B, 0)

    backend = get_backend(A)
    kernel! = tr_kernel!(backend)

    kernel!(B, A, ndrange = size(A)[1:2])
    KernelAbstractions.synchronize(backend)

    return B
end

end
