#=
Functions which manipulates QuantumObject
=#

export ket2dm
export expect, variance
export sparse_to_dense, dense_to_sparse
export vec2mat, mat2vec

@doc raw"""
    ket2dm(ψ::QuantumObject)

Transform the ket state ``\ket{\psi}`` into a pure density matrix ``\hat{\rho} = \dyad{\psi}``.
"""
ket2dm(ψ::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} = ψ * ψ'

ket2dm(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = ρ

@doc raw"""
    expect(O::AbstractQuantumObject, ψ::Union{QuantumObject,Vector{QuantumObject}})

Expectation value of the [`Operator`](@ref) `O` with the state `ψ`. The state can be a [`Ket`](@ref), [`Bra`](@ref) or [`Operator`](@ref).

If `ψ` is a [`Ket`](@ref) or [`Bra`](@ref), the function calculates ``\langle\psi|\hat{O}|\psi\rangle``.

If `ψ` is a density matrix ([`Operator`](@ref)), the function calculates ``\textrm{Tr} \left[ \hat{O} \hat{\psi} \right]``

The function returns a real number if `O` is of `Hermitian` type or `Symmetric` type, and returns a complex number otherwise. You can make an operator `O` hermitian by using `Hermitian(O)`.

Note that `ψ` can also be given as a list of [`QuantumObject`](@ref), it returns a list of expectation values.

# Examples

```
julia> ψ = 1 / √2 * (fock(10,2) + fock(10,4));

julia> a = destroy(10);

julia> expect(a' * a, ψ) |> round
3.0 + 0.0im

julia> expect(Hermitian(a' * a), ψ) |> round
3.0
```
"""
function expect(
    O::AbstractQuantumObject{DT1,OperatorQuantumObject},
    ψ::QuantumObject{DT2,KetQuantumObject},
) where {DT1,DT2}
    return dot(ψ.data, O.data, ψ.data)
end
function expect(
    O::AbstractQuantumObject{DT1,OperatorQuantumObject},
    ψ::QuantumObject{DT2,BraQuantumObject},
) where {DT1,DT2}
    return expect(O, ψ')
end
function expect(
    O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    ρ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    return tr(O * ρ)
end
function expect(
    O::QuantumObject{<:Union{<:Hermitian{TF},<:Symmetric{TR}},OperatorQuantumObject},
    ψ::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {TF<:Number,TR<:Real,T2}
    return real(dot(ψ.data, O.data, ψ.data))
end
function expect(
    O::QuantumObject{<:Union{<:Hermitian{TF},<:Symmetric{TR}},OperatorQuantumObject},
    ψ::QuantumObject{<:AbstractArray{T2},BraQuantumObject},
) where {TF<:Number,TR<:Real,T2}
    return real(expect(O, ψ'))
end
function expect(
    O::QuantumObject{<:Union{<:Hermitian{TF},<:Symmetric{TR}},OperatorQuantumObject},
    ρ::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {TF<:Number,TR<:Real,T2}
    return real(tr(O * ρ))
end
function expect(O::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, ρ::Vector{<:QuantumObject}) where {T1}
    _expect = _ρ -> expect(O, _ρ)
    return _expect.(ρ)
end

@doc raw"""
    variance(O::QuantumObject, ψ::Union{QuantumObject,Vector{QuantumObject}})

Variance of the [`Operator`](@ref) `O`: ``\langle\hat{O}^2\rangle - \langle\hat{O}\rangle^2``,

where ``\langle\hat{O}\rangle`` is the expectation value of `O` with the state `ψ` (see also [`expect`](@ref)), and the state `ψ` can be a [`Ket`](@ref), [`Bra`](@ref) or [`Operator`](@ref).

The function returns a real number if `O` is hermitian, and returns a complex number otherwise.

Note that `ψ` can also be given as a list of [`QuantumObject`](@ref), it returns a list of expectation values.
"""
variance(O::QuantumObject{DT1,OperatorQuantumObject}, ψ::QuantumObject{DT2}) where {DT1,DT2} =
    expect(O^2, ψ) - expect(O, ψ)^2
variance(O::QuantumObject{DT1,OperatorQuantumObject}, ψ::Vector{<:QuantumObject}) where {DT1} =
    expect(O^2, ψ) .- expect(O, ψ) .^ 2

@doc raw"""
    sparse_to_dense(A::QuantumObject)

Converts a sparse QuantumObject to a dense QuantumObject.
"""
sparse_to_dense(A::QuantumObject{<:AbstractVecOrMat}) = QuantumObject(sparse_to_dense(A.data), A.type, A.dims)
sparse_to_dense(A::MT) where {MT<:AbstractSparseArray} = Array(A)
for op in (:Transpose, :Adjoint)
    @eval sparse_to_dense(A::$op{T,<:AbstractSparseMatrix}) where {T<:BlasFloat} = Array(A)
end
sparse_to_dense(A::MT) where {MT<:AbstractArray} = A

sparse_to_dense(::Type{T}, A::AbstractSparseArray) where {T<:Number} = Array{T}(A)
sparse_to_dense(::Type{T1}, A::AbstractArray{T2}) where {T1<:Number,T2<:Number} = Array{T1}(A)
sparse_to_dense(::Type{T}, A::AbstractArray{T}) where {T<:Number} = A

function sparse_to_dense(::Type{M}) where {M<:SparseMatrixCSC}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return Matrix{par[1]}
end

sparse_to_dense(::Type{M}) where {M<:AbstractMatrix} = M

@doc raw"""
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

@doc raw"""
    kron(A::AbstractQuantumObject, B::AbstractQuantumObject, ...)

Returns the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) ``\hat{A} \otimes \hat{B} \otimes \cdots``.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> kron(a, a)
Quantum Object:   type=Operator   dims=[20, 20]   size=(400, 400)   ishermitian=false
400×400 SparseMatrixCSC{ComplexF64, Int64} with 361 stored entries:
⠀⠀⠘⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠦
```
"""
function LinearAlgebra.kron(
    A::AbstractQuantumObject{DT1,OpType},
    B::AbstractQuantumObject{DT2,OpType},
) where {DT1,DT2,OpType<:Union{KetQuantumObject,BraQuantumObject,OperatorQuantumObject}}
    if A isa QuantumObjectEvolution || B isa QuantumObjectEvolution
        return QuantumObjectEvolution(kron(A.data, B.data), A.type, vcat(A.dims, B.dims))
    end
    return QuantumObject(kron(A.data, B.data), A.type, vcat(A.dims, B.dims))
end
LinearAlgebra.kron(A::AbstractQuantumObject) = A
function LinearAlgebra.kron(A::Vector{<:AbstractQuantumObject})
    @warn "`tensor(A)` or `kron(A)` with `A` is a `Vector` can hurt performance. Try to use `tensor(A...)` or `kron(A...)` instead."
    return kron(A...)
end

@doc raw"""
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

@doc raw"""
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
