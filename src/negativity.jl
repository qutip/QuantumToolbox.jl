@doc raw"""
    partial_transpose(ρ::QuantumObject, mask::Vector{Bool})

Return the partial transpose of a density matrix ``\rho``, where `mask` is an array/vector with length that equals the length of `ρ.dims`. The elements in `mask` are boolean (`true` or `false`) which indicates whether or not the corresponding subsystem should be transposed.
"""
function partial_transpose(ρ::QuantumObject{::T, OperatorQuantumObject}, mask::Vector{Bool}) where T
    if length(mask) != length(ρ.dims) 
        error("The length of \`mask\` should be equal to the length of \`ρ.dims\`.")
    end
    return _partial_transpose(ρ, mask)
end

function _partial_transpose(ρ::QuantumObject{<:AbstractMatrix, OperatorQuantumObject}, mask::Vector{Bool})
    mask2 = [1 + Int(i) for i in mask]
    # mask2 has elements with values equal to 1 or 2
    # 1 - the subsystem don't need to be transposed
    # 2 - the subsystem need be transposed

    nsys = length(mask2)
    pt_dims = reshape(Vector(1:(2 * nsys)), (nsys, 2))
    pt_idx  = [
        [pt_dims[n,     mask2[n]] for n in 1:nsys]; # origin   value in mask2
        [pt_dims[n, 3 - mask2[n]] for n in 1:nsys]  # opposite value in mask2 (1 -> 2, and 2 -> 1)
    ]
    return QuantumObject(
        reshape(permutedims(reshape(ρ.data, (ρ.dims..., ρ.dims...)), pt_idx), size(ρ)),
        OperatorQuantumObject,
        ρ.dims
    )
end