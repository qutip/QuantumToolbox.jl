#=
Functions which manipulates QuantumObject
=#

export ket2dm
export expect
export sparse_to_dense, dense_to_sparse
export vec2mat, mat2vec

@doc raw"""
    ket2dm(ψ::QuantumObject)

Transform the ket state ``\ket{\psi}`` into a pure density matrix ``\hat{\rho} = \dyad{\psi}``.
"""
ket2dm(ψ::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} = ψ * ψ'

ket2dm(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = ρ

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
function expect(
    O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2}
    ψd = ψ.data
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od, ψd)) : dot(ψd, Od, ψd)
end
function expect(
    O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ψ::QuantumObject{<:AbstractArray{T2},BraQuantumObject},
) where {T1,T2}
    ψd = ψ.data'
    Od = O.data
    return ishermitian(O) ? real(dot(ψd, Od, ψd)) : dot(ψd, Od, ψd)
end
function expect(
    O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ρ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    return ishermitian(O) ? real(tr(O * ρ)) : tr(O * ρ)
end

"""
    sparse_to_dense(A::QuantumObject)

Converts a sparse QuantumObject to a dense QuantumObject.
"""
sparse_to_dense(A::QuantumObject{<:AbstractVecOrMat}) = QuantumObject(sparse_to_dense(A.data), A.type, A.dims)
sparse_to_dense(A::MT) where {MT<:AbstractSparseMatrix} = Array(A)
for op in (:Transpose, :Adjoint)
    @eval sparse_to_dense(A::$op{T,<:AbstractSparseMatrix}) where {T<:BlasFloat} = Array(A)
end
sparse_to_dense(A::MT) where {MT<:AbstractArray} = A

function sparse_to_dense(::Type{M}) where {M<:SparseMatrixCSC}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return Matrix{par[1]}
end

sparse_to_dense(::Type{M}) where {M<:AbstractMatrix} = M

"""
    dense_to_sparse(A::QuantumObject)

Converts a dense QuantumObject to a sparse QuantumObject.
"""
dense_to_sparse(A::QuantumObject{<:AbstractVecOrMat}, tol::Real = 1e-10) =
    QuantumObject(dense_to_sparse(A.data, tol), A.type, A.dims)
function dense_to_sparse(A::MT, tol::Real = 1e-10) where {MT<:AbstractMatrix}
    idxs = findall(@. abs(A) > tol)
    row_indices = getindex.(idxs, 1)
    col_indices = getindex.(idxs, 2)
    vals = getindex(A, idxs)
    return sparse(row_indices, col_indices, vals, size(A)...)
end
function dense_to_sparse(A::VT, tol::Real = 1e-10) where {VT<:AbstractVector}
    idxs = findall(@. abs(A) > tol)
    vals = getindex(A, idxs)
    return sparsevec(idxs, vals, length(A))
end

"""
    vec2mat(A::AbstractVector)

Converts a vector to a matrix.
"""
function vec2mat(A::AbstractVector)
    newsize = isqrt(length(A))
    return reshape(A, newsize, newsize)
end

@doc raw"""
    vec2mat(A::QuantumObject)

Convert a quantum object from vector ([`OperatorKetQuantumObject`](@ref)-type) to matrix ([`OperatorQuantumObject`](@ref)-type)
"""
vec2mat(A::QuantumObject{<:AbstractArray{T},OperatorKetQuantumObject}) where {T} =
    QuantumObject(vec2mat(A.data), Operator, A.dims)

@doc raw"""
    mat2vec(A::QuantumObject)

Convert a quantum object from matrix ([`OperatorQuantumObject`](@ref)-type) to vector ([`OperatorKetQuantumObject`](@ref)-type)
"""
mat2vec(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} =
    QuantumObject(mat2vec(A.data), OperatorKet, A.dims)

"""
    mat2vec(A::AbstractMatrix)

Converts a matrix to a vector.
"""
mat2vec(A::MT) where {MT<:AbstractMatrix} = vec(A) # reshape(A, :)
function mat2vec(A::MT) where {MT<:AbstractSparseMatrix}
    i, j, v = findnz(A)
    return sparsevec(i .+ (j .- 1) .* size(A, 1), v, prod(size(A)))
end
for op in (:Transpose, :Adjoint)
    @eval mat2vec(A::$op{T,<:AbstractSparseMatrix}) where {T<:BlasFloat} = mat2vec(sparse(A))
    @eval mat2vec(A::$op{T,<:AbstractMatrix}) where {T<:BlasFloat} = mat2vec(Matrix(A))
end

function mat2vec(::Type{M}) where {M<:DenseMatrix}
    T = hasproperty(M, :body) ? M.body : M
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $M is not supported.")
    if npar == 2
        S = T.name.wrapper{par[1],1}
    else
        S = T.name.wrapper{par[1],1,par[3]}
    end
    return S
end

function mat2vec(::Type{M}) where {M<:SparseMatrixCSC}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return SparseVector{par[1],par[2]}
end

function mat2vec(
    ::Type{M},
) where {M<:Union{Adjoint{<:BlasFloat,<:SparseMatrixCSC},Transpose{<:BlasFloat,<:SparseMatrixCSC}}}
    T = M.parameters[2]
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return SparseVector{par[1],par[2]}
end
