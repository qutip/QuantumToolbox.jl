"""
    row_major_reshape(Q::AbstractArray, shapes...)

Reshapes `Q` in the row-major order, as numpy. 
"""
row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} = PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner=length(y)), length(y), length(x))
    Y = repeat(y, outer=(1, length(x)))
    X, Y
end

"""
    sparse_to_dense(A::QuantumObject)

Converts a sparse QuantumObject to a dense QuantumObject.
"""
sparse_to_dense(A::QuantumObject{<:AbstractArray{T}}) where T = QuantumObject(sparse_to_dense(A.data), A.type, A.dims)
sparse_to_dense(A::AbstractSparseArray) = Array(A)
sparse_to_dense(A::AbstractArray) = A

"""
    dense_to_sparse(A::QuantumObject)

Converts a dense QuantumObject to a sparse QuantumObject.
"""
dense_to_sparse(A::QuantumObject{<:AbstractArray{T}}) where T = QuantumObject(dense_to_sparse(A.data), A.type, A.dims)
function dense_to_sparse(A::AbstractMatrix, tol::Real=1e-10)
    idxs = findall(abs.(A) .> tol)
    row_indices = getindex.(idxs, 1)
    col_indices = getindex.(idxs, 2)
    vals = getindex(A, idxs)
    return sparse(row_indices, col_indices, vals, size(A)...)
end
function dense_to_sparse(A::Vector, tol::Real=1e-10)
    idxs = findall(abs.(A) .> tol)
    return sparse(idxs, vals, length(A))
end

@doc raw"""
    gaussian(x, μ::Real, σ::Real)

Returns the gaussian function ``\exp \left[- \frac{(x - \mu)^2}{2 \sigma^2} \right]``,
where ``\mu`` and ``\sigma^2`` are the mean and the variance respectively.
"""
gaussian(x::AbstractVector{T}, μ::Real, σ::Real) where {T} = @. exp(-0.5 * (x - μ)^2 / σ^2)

gaussian(x::Real, μ::Real, σ::Real) = exp(-0.5 * (x - μ)^2 / σ^2)

@doc raw"""
    ptrace(QO::QuantumObject, sel::Vector{Int})

[Partial trace](https://en.wikipedia.org/wiki/Partial_trace) of a quantum state `QO` leaving only the dimensions
with the indices present in the `sel` vector.

# Examples
Two qubits in the state ``\ket{\psi} = \ket{e,g}``:
```
julia> ψ = kron(fock(2,0), fock(2,1))
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> ptrace(ψ, [2])
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im
```

or in an entangled state ``\ket{\psi} = \frac{1}{\sqrt{2}} \left( \ket{e,e} + \ket{g,g} \right)``:
```
julia> ψ = 1 / √2 * (kron(fock(2,0), fock(2,0)) + kron(fock(2,1), fock(2,1)))
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> ptrace(ψ, [1])
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im
```
"""
function ptrace(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::Vector{Int}) where
{T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    rd = QO.dims
    nd = length(rd)
    dkeep = rd[sel]
    qtrace = filter(e -> e ∉ sel, 1:nd)
    dtrace = rd[qtrace]

    nd == 1 && return QO

    if isket(QO) || isbra(QO)
        vmat = row_major_reshape(QO.data, prod(rd), 1)
        vmat = row_major_reshape(vmat, rd...)
        topermute = Int64[]
        append!(topermute, sel)
        append!(topermute, qtrace)
        reverse!(topermute)
        vmat = PermutedDimsArray(vmat, topermute)
        vmat = row_major_reshape(vmat, prod(dkeep), prod(dtrace))
        return QuantumObject(vmat * vmat', OperatorQuantumObject, dkeep)
    elseif isoper(QO)
        ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
        topermute = Int64[]
        append!(topermute, qtrace)
        append!(topermute, [nd + q for q in qtrace])
        append!(topermute, sel)
        append!(topermute, [nd + q for q in sel])
        reverse!(topermute)
        ρmat = PermutedDimsArray(ρmat, topermute)
        ρmat = row_major_reshape(ρmat, prod(dtrace), prod(dtrace), prod(dkeep), prod(dkeep))
        dims = size(ρmat)
        res = [tr(@view(ρmat[:, :, i, j])) for i in 1:dims[3] for j in 1:dims[4]]
        return QuantumObject(row_major_reshape(res, dims[3:length(dims)]...), OperatorQuantumObject, dkeep)
    end
end

@doc raw"""
    entropy_vn(ρ::QuantumObject; base::Int=0, tol::Real=1e-15)

Calculates the [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy)
``S = - \Tr \left[ \hat{\rho} \log \left( \hat{\rho} \right) \right]`` where ``\hat{\rho}``
is the density matrix of the system.

The `base` parameter specifies the base of the logarithm to use, and when using the default value 0,
the natural logarithm is used. The `tol` parameter
describes the absolute tolerance for detecting the zero-valued eigenvalues of the density
matrix ``\hat{\rho}``.

# Examples

Pure state:
```
julia> ψ = fock(2,0)
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im

julia> ρ = ket2dm(ψ)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> entropy_vn(ρ, base=2)
-0.0
```

Mixed state:
```
julia> ρ = 1 / 2 * ( ket2dm(fock(2,0)) + ket2dm(fock(2,1)) )
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im

julia> entropy_vn(ρ, base=2)
1.0
```
"""
function entropy_vn(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}; base::Int=0, tol::Real=1e-15) where {T}
    vals = eigvals(ρ)
    indexes = abs.(vals) .> tol
    1 ∉ indexes && return 0 
    nzvals = vals[indexes]
    logvals = base != 0 ? log.(base, Complex.(nzvals)) : log.(Complex.(nzvals))
    return -real(sum(nzvals .* logvals))
end

"""
    entanglement(QO::QuantumObject, sel::Vector)

Calculates the entanglement by doing the partial trace of `QO`, selecting only the dimensions
with the indices contained in the `sel` vector, and then using the Von Neumann entropy [`entropy_vn`](@ref). 
"""
function entanglement(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::Vector{Int}) where
{T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}

    ψ = normalize(QO)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end

@doc raw"""
    expect(O::QuantumObject, ψ::QuantumObject)

Expectation value of the operator `O` with the state `ψ`. The latter
can be a [`KetQuantumObject`](@ref), [`BraQuantumObject`](@ref) or a [`OperatorQuantumObject`](@ref).
If `ψ` is a density matrix, the function calculates ``\Tr \left[ \hat{O} \hat{\psi} \right]``, while if `ψ`
is a state, the function calculates ``\mel{\psi}{\hat{O}}{\psi}``.

The function returns a real number if the operator is hermitian, and returns a complex number otherwise.

# Examples

```
julia> ψ = 1 / √2 * (fock(10,2) + fock(10,4));

julia> a = destroy(10);

julia> expect(a' * a, ψ) ≈ 3
true
```
"""
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T2},KetQuantumObject}) where {T1,T2}
    ψd = ψ.data
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od * ψd)) : dot(ψd, Od * ψd)
end
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ψ::QuantumObject{<:AbstractArray{T2},BraQuantumObject}) where {T1,T2}
    ψd = ψ.data'
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od * ψd)) : dot(ψd, Od * ψd)
end
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ρ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject}) where {T1,T2}
    return ishermitian(O) ? real(tr(O * ρ)) : tr(O * ρ)
end

@doc raw"""
    get_coherence(ψ::QuantumObject)

Get the coherence value ``\alpha`` by measuring the expectation value of the destruction
operator ``\hat{a}`` on the state ``\ket{\psi}``.

It returns both ``\alpha`` and the state 
``\ket{\delta_\psi} = \exp ( \bar{\alpha} \hat{a} - \alpha \hat{a}^\dagger )``. The
latter corresponds to the quantum fulctuations around the coherent state ``\ket{\alpha}``.
"""
function get_coherence(ψ::QuantumObject{<:AbstractArray{T}, StateOpType}) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    a = destroy(size(ψ,1))
    α = expect(a, ψ)
    D = exp(α*a' - conj(α)*a)

    α, D' * ψ
end

@doc raw"""
    n_th(ω::Number, T::Real)

Gives the mean number of excitations in a mode with frequency ω at temperature T:
``n_{\rm th} (\omega, T) = \frac{1}{e^{\omega/T} - 1}``
"""
function n_th(ω::Real, T::Real)::Float64
    (T == 0 || ω == 0) && return 0.0
    abs(ω / T) > 50 && return 0.0
    return 1 / (exp(ω / T) - 1)
end
