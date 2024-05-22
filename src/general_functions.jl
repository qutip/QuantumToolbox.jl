export get_data, get_coherence, expect, ptrace
export mat2vec, vec2mat
export entropy_vn, entanglement, tracedist
export gaussian, n_th

export row_major_reshape, tidyup, tidyup!, meshgrid, sparse_to_dense, dense_to_sparse

"""
    row_major_reshape(Q::AbstractArray, shapes...)

Reshapes `Q` in the row-major order, as numpy.
"""
row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} =
    PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    return X, Y
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
    tidyup(A::QuantumObject, tol::Real=1e-14)

Removes those elements of a QuantumObject `A` whose absolute value is less than `tol`.
"""
tidyup(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} =
    QuantumObject(tidyup(A.data, tol), A.type, A.dims)
tidyup(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} = @. T(abs(A) > tol) * A
tidyup(A::AbstractSparseMatrix{T}, tol::T2 = 1e-14) where {T,T2<:Real} = droptol!(copy(A), tol)

"""
    tidyup!(A::QuantumObject, tol::Real=1e-14)

Removes those elements of a QuantumObject `A` whose absolute value is less than `tol`.
In-place version of [`tidyup`](#tidyup).
"""
tidyup!(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} = (tidyup!(A.data, tol); A)
tidyup!(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} = @. A = T(abs(A) > tol) * A
tidyup!(A::AbstractSparseMatrix{T}, tol::T2 = 1e-14) where {T,T2<:Real} = droptol!(A, tol)

"""
    get_data(A::QuantumObject)

Returns the data of a QuantumObject.
"""
get_data(A::QuantumObject) = A.data

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
    gaussian(x::Number, μ::Number, σ::Number)

Returns the gaussian function ``\exp \left[- \frac{(x - \mu)^2}{2 \sigma^2} \right]``,
where ``\mu`` and ``\sigma^2`` are the mean and the variance respectively.
"""
gaussian(x::Number, μ::Number, σ::Number) = exp(-(x - μ)^2 / (2 * σ^2))

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
function ptrace(QO::QuantumObject{<:AbstractArray{T1},KetQuantumObject}, sel::Vector{T2}) where {T1,T2<:Integer}
    length(QO.dims) == 1 && return QO

    ρtr, dkeep = _ptrace_ket(QO.data, QO.dims, sel)
    return QuantumObject(ρtr, dims = dkeep)
end

ptrace(QO::QuantumObject{<:AbstractArray{T1},BraQuantumObject}, sel::Vector{T2}) where {T1,T2<:Integer} =
    ptrace(QO', sel)

function ptrace(QO::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject}, sel::Vector{T2}) where {T1,T2<:Integer}
    length(QO.dims) == 1 && return QO

    ρtr, dkeep = _ptrace_oper(QO.data, QO.dims, sel)
    return QuantumObject(ρtr, dims = dkeep)
end
ptrace(QO::QuantumObject, sel::Int) = ptrace(QO, [sel])

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
function entropy_vn(
    ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject};
    base::Int = 0,
    tol::Real = 1e-15,
) where {T}
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
function entanglement(
    QO::QuantumObject{<:AbstractArray{T},OpType},
    sel::Vector{Int},
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    ψ = normalize(QO)
    ρ_tr = ptrace(ψ, sel)
    entropy = entropy_vn(ρ_tr)
    return (entropy > 0) * entropy
end
entanglement(QO::QuantumObject, sel::Int) = entanglement(QO, [sel])

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

@doc raw"""
    get_coherence(ψ::QuantumObject)

Get the coherence value ``\alpha`` by measuring the expectation value of the destruction
operator ``\hat{a}`` on the state ``\ket{\psi}``.

It returns both ``\alpha`` and the state 
``\ket{\delta_\psi} = \exp ( \bar{\alpha} \hat{a} - \alpha \hat{a}^\dagger )``. The
latter corresponds to the quantum fulctuations around the coherent state ``\ket{\alpha}``.
"""
function get_coherence(
    ψ::QuantumObject{<:AbstractArray{T},StateOpType},
) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    a = destroy(size(ψ, 1))
    α = expect(a, ψ)
    D = exp(α * a' - conj(α) * a)

    return α, D' * ψ
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

@doc raw"""
    tracedist(ρ::QuantumObject, σ::QuantumObject)

Calculates the [trace distance](https://en.wikipedia.org/wiki/Trace_distance) between two [`QuantumObject`](@ref):
``T(\rho, \sigma) = \frac{1}{2} \lVert \rho - \sigma \rVert_1``

Note that `ρ` and `σ` must be either [`Ket`](@ref) or [`Operator`](@ref).
"""
tracedist(
    ρ::QuantumObject{<:AbstractArray{T1},ObjType1},
    σ::QuantumObject{<:AbstractArray{T2},ObjType2},
) where {
    T1,
    T2,
    ObjType1<:Union{KetQuantumObject,OperatorQuantumObject},
    ObjType2<:Union{KetQuantumObject,OperatorQuantumObject},
} = norm(ket2dm(ρ) - ket2dm(σ), 1) / 2

function _ptrace_ket(QO::AbstractArray{T1}, dims::Vector{<:Integer}, sel::Vector{T2}) where {T1,T2<:Integer}
    rd = dims
    nd = length(rd)

    nd == 1 && return QO, rd

    dkeep = rd[sel]
    qtrace = filter!(e -> e ∉ sel, Vector(1:nd))
    dtrace = @view(rd[qtrace])

    vmat = reshape(QO, reverse(rd)...)
    topermute = nd + 1 .- vcat(sel, qtrace)
    vmat = PermutedDimsArray(vmat, topermute)
    vmat = reshape(vmat, prod(dkeep), prod(dtrace))

    return vmat * vmat', dkeep
end

function _ptrace_oper(QO::AbstractArray{T1}, dims::Vector{<:Integer}, sel::Vector{T2}) where {T1,T2<:Integer}
    rd = dims
    nd = length(rd)

    nd == 1 && return QO, rd

    dkeep = rd[sel]
    qtrace = filter!(e -> e ∉ sel, Vector(1:nd))
    dtrace = @view(rd[qtrace])

    ρmat = reshape(QO, reverse!(repeat(rd, 2))...)
    topermute = 2 * nd + 1 .- vcat(qtrace, qtrace .+ nd, sel, sel .+ nd)
    reverse!(topermute)
    ρmat = PermutedDimsArray(ρmat, topermute)

    ## TODO: Check if it works always

    # ρmat = row_major_reshape(ρmat, prod(dtrace), prod(dtrace), prod(dkeep), prod(dkeep))
    # res = dropdims(mapslices(tr, ρmat, dims=(1,2)), dims=(1,2))
    ρmat = reshape(ρmat, prod(dkeep), prod(dkeep), prod(dtrace), prod(dtrace))
    res = dropdims(mapslices(tr, ρmat, dims = (3, 4)), dims = (3, 4))

    return res, dkeep
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

@doc raw"""
    mat2vec(A::QuantumObject)

Convert a quantum object from matrix ([`OperatorQuantumObject`](@ref)-type) to vector ([`OperatorKetQuantumObject`](@ref)-type)
"""
mat2vec(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} =
    QuantumObject(mat2vec(A.data), OperatorKet, A.dims)
