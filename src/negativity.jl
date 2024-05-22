export negativity, partial_transpose

@doc raw"""
    negativity(ρ::QuantumObject, subsys::Int; logarithmic::Bool=false)

Compute the [negativity](https://en.wikipedia.org/wiki/Negativity_(quantum_mechanics)) ``N(\rho) = \frac{\Vert \rho^{\Gamma}\Vert_1 - 1}{2}``  
where ``\rho^{\Gamma}`` is the partial transpose of ``\rho`` with respect to the subsystem,  
and ``\Vert X \Vert_1=\textrm{Tr}\sqrt{X^\dagger X}`` is the trace norm.

# Arguments
- `ρ::QuantumObject`: The density matrix (`ρ.type` must be [`OperatorQuantumObject`](@ref)).
- `subsys::Int`: an index that indicates which subsystem to compute the negativity for.
- `logarithmic::Bool`: choose whether to calculate logarithmic negativity or not. Default as `false`

# Returns
- `N::Real`: The value of negativity.

# Examples

```
julia> Ψ = 1 / √2 * ( basis(2,0) ⊗ basis(2,0) + basis(2,1) ⊗ basis(2,1) )
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> ρ = ket2dm(Ψ)
Quantum Object:   type=Operator   dims=[2, 2]   size=(4, 4)   ishermitian=true
4×4 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im

julia> negativity(ρ, 2)
0.4999999999999998
```
"""
function negativity(ρ::QuantumObject, subsys::Int; logarithmic::Bool = false)
    mask = fill(false, length(ρ.dims))
    try
        mask[subsys] = true
    catch
        error("Invalid index of subsys: $subsys")
    end

    ρ_pt = partial_transpose(ρ, mask)
    tr_norm = norm(ρ_pt, 1) # trace norm

    if logarithmic
        return log2(tr_norm)
    else
        return (tr_norm - 1) / 2
    end
end

@doc raw"""
    partial_transpose(ρ::QuantumObject, mask::Vector{Bool})

Return the partial transpose of a density matrix ``\rho``, where `mask` is an array/vector with length that equals the length of `ρ.dims`. The elements in `mask` are boolean (`true` or `false`) which indicates whether or not the corresponding subsystem should be transposed.

# Arguments
- `ρ::QuantumObject`: The density matrix (`ρ.type` must be [`OperatorQuantumObject`](@ref)).
- `mask::Vector{Bool}`: A boolean vector selects which subsystems should be transposed.

# Returns
- `ρ_pt::QuantumObject`: The density matrix with the selected subsystems transposed.
"""
function partial_transpose(ρ::QuantumObject{T,OperatorQuantumObject}, mask::Vector{Bool}) where {T}
    if length(mask) != length(ρ.dims)
        error("The length of \`mask\` should be equal to the length of \`ρ.dims\`.")
    end
    return _partial_transpose(ρ, mask)
end

# for dense matrices
function _partial_transpose(ρ::QuantumObject{<:AbstractArray,OperatorQuantumObject}, mask::Vector{Bool})
    mask2 = [1 + Int(i) for i in mask]
    # mask2 has elements with values equal to 1 or 2
    #   1 - the subsystem don't need to be transposed
    #   2 - the subsystem need be transposed

    nsys = length(mask2)
    pt_dims = reshape(Vector(1:(2*nsys)), (nsys, 2))
    pt_idx = [
        [pt_dims[n, mask2[n]] for n in 1:nsys] # origin   value in mask2
        [pt_dims[n, 3-mask2[n]] for n in 1:nsys]  # opposite value in mask2 (1 -> 2, and 2 -> 1)
    ]
    return QuantumObject(
        reshape(permutedims(reshape(ρ.data, (ρ.dims..., ρ.dims...)), pt_idx), size(ρ)),
        Operator,
        ρ.dims,
    )
end

# for sparse matrices
function _partial_transpose(ρ::QuantumObject{<:AbstractSparseArray,OperatorQuantumObject}, mask::Vector{Bool})
    M, N = size(ρ)
    dimsTuple = Tuple(ρ.dims)
    colptr = ρ.data.colptr
    rowval = ρ.data.rowval
    nzval = ρ.data.nzval
    len = length(nzval)

    # for partial transposed data
    I_pt = Vector{Int}(undef, len)
    J_pt = Vector{Int}(undef, len)
    V_pt = Vector{eltype(ρ)}(undef, len)

    n = 0
    for j in 1:(length(colptr)-1)
        for p in colptr[j]:(colptr[j+1]-1)
            n += 1
            i = rowval[p]
            if i == j
                I_pt[n] = i
                J_pt[n] = j
            else
                ket_pt = [Base._ind2sub(dimsTuple, i)...]
                bra_pt = [Base._ind2sub(dimsTuple, j)...]
                for sys in findall(m -> m, mask)
                    @inbounds ket_pt[sys], bra_pt[sys] = bra_pt[sys], ket_pt[sys]
                end
                I_pt[n] = Base._sub2ind(dimsTuple, ket_pt...)
                J_pt[n] = Base._sub2ind(dimsTuple, bra_pt...)
            end
            V_pt[n] = nzval[p]
        end
    end

    return QuantumObject(sparse(I_pt, J_pt, V_pt, M, N), Operator, ρ.dims)
end
