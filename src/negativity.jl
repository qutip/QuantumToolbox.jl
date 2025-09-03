export negativity, partial_transpose

@doc raw"""
    negativity(ρ::QuantumObject, subsys::Int; logarithmic::Bool=false)

Compute the [negativity](https://en.wikipedia.org/wiki/Negativity_(quantum_mechanics)) ``N(\hat{\rho}) = \frac{\Vert \hat{\rho}^{\Gamma}\Vert_1 - 1}{2}``  
where ``\hat{\rho}^{\Gamma}`` is the partial transpose of ``\hat{\rho}`` with respect to the subsystem,  
and ``\Vert \hat{X} \Vert_1=\textrm{Tr}\sqrt{\hat{X}^\dagger \hat{X}}`` is the trace norm.

# Arguments
- `ρ::QuantumObject`: The density matrix (`ρ.type` must be [`Operator`](@ref)).
- `subsys::Int`: an index that indicates which subsystem to compute the negativity for.
- `logarithmic::Bool`: choose whether to calculate logarithmic negativity or not. Default as `false`

# Returns
- `N::Real`: The value of negativity.

# Examples

```jldoctest
julia> Ψ = bell_state(0, 0)

Quantum Object:   type=Ket()   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> ρ = ket2dm(Ψ)

Quantum Object:   type=Operator()   dims=[2, 2]   size=(4, 4)   ishermitian=true
4×4 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im

julia> round(negativity(ρ, 2), digits=2)
0.5
```
"""
function negativity(ρ::QuantumObject, subsys::Int; logarithmic::Bool = false)
    mask = fill(false, length(ρ.dimensions))
    try
        mask[subsys] = true
    catch
        throw(ArgumentError("Invalid index of subsys: $subsys"))
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
- `ρ::QuantumObject`: The density matrix (`ρ.type` must be [`Operator`](@ref)).
- `mask::Vector{Bool}`: A boolean vector selects which subsystems should be transposed.

# Returns
- `ρ_pt::QuantumObject`: The density matrix with the selected subsystems transposed.
"""
function partial_transpose(ρ::QuantumObject{Operator}, mask::Vector{Bool})
    any(s -> s isa EnrSpace, ρ.dimensions.to) && throw(ArgumentError("partial_transpose does not support EnrSpace"))

    (length(mask) != length(ρ.dimensions)) &&
        throw(ArgumentError("The length of \`mask\` should be equal to the length of \`ρ.dims\`."))

    isa(ρ.dimensions, GeneralDimensions) &&
        (get_dimensions_to(ρ) != get_dimensions_from(ρ)) &&
        throw(ArgumentError("Invalid partial transpose for dims = $(_get_dims_string(ρ.dimensions))"))

    return _partial_transpose(ρ, mask)
end

# for dense matrices
function _partial_transpose(ρ::QuantumObject{Operator}, mask::Vector{Bool})
    nsys = length(mask)
    mask2 = reverse([mask[s] ? 2 : 1 for s in 1:nsys])
    # mask2 has elements with values equal to 1 or 2
    #   1 - the subsystem (in reversed order) don't need to be transposed
    #   2 - the subsystem (in reversed order) need to be transposed

    dims_rev = reverse(dimensions_to_dims(get_dimensions_to(ρ)))
    pt_dims = reshape(Vector(1:(2*nsys)), (nsys, 2))
    pt_idx = [
        [pt_dims[n, mask2[n]] for n in 1:nsys]   # origin   value in mask2
        [pt_dims[n, 3-mask2[n]] for n in 1:nsys] # opposite value in mask2 (1 -> 2, and 2 -> 1)
    ]
    return QuantumObject(
        reshape(permutedims(reshape(ρ.data, (dims_rev..., dims_rev...)), pt_idx), size(ρ)),
        Operator(),
        Dimensions(ρ.dimensions.to),
    )
end

# for sparse matrices
function _partial_transpose(
    ρ::QuantumObject{Operator,DimsType,<:AbstractSparseArray},
    mask::Vector{Bool},
) where {DimsType<:AbstractDimensions}
    M, N = size(ρ)
    dims_rev = reverse(Tuple(dimensions_to_dims(get_dimensions_to(ρ))))
    mask_rev = reverse(mask)
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
                ket_pt = [Base._ind2sub(dims_rev, i)...]
                bra_pt = [Base._ind2sub(dims_rev, j)...]
                for sys in findall(m -> m, mask_rev)
                    @inbounds ket_pt[sys], bra_pt[sys] = bra_pt[sys], ket_pt[sys]
                end
                I_pt[n] = Base._sub2ind(dims_rev, ket_pt...)
                J_pt[n] = Base._sub2ind(dims_rev, bra_pt...)
            end
            V_pt[n] = nzval[p]
        end
    end

    return QuantumObject(sparse(I_pt, J_pt, V_pt, M, N), Operator(), ρ.dimensions)
end
