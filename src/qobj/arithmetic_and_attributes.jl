#=
Arithmetic and Attributes for QuantumObject
    - extend most of the useful functions in LinearAlgebra for QuantumObject
    - export most of the attribute functions in "Python Qobj class"
=#

export trans, dag, dagger, matrix_element
export sqrtm, sinm, cosm
export ptrace
export tidyup, tidyup!
export get_data, get_coherence

#    Broadcasting
Base.broadcastable(x::QuantumObject) = x.data
for op in (:(+), :(-), :(*), :(/), :(^))
    @eval begin
        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::QuantumObject)
            return QuantumObject(broadcast($op, x.data, y.data), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::Number)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::Number, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::AbstractArray)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::AbstractArray, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dims)
        end
    end
end

for op in (:(+), :(-), :(*))
    @eval begin
        function LinearAlgebra.$op(
            A::QuantumObject{<:AbstractArray{T1},OpType},
            B::QuantumObject{<:AbstractArray{T2},OpType},
        ) where {T1,T2,OpType<:QuantumObjectType}
            A.dims != B.dims &&
                throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
            return QuantumObject($(op)(A.data, B.data), A.type, A.dims)
        end
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T}}) where {T} = QuantumObject($(op)(A.data), A.type, A.dims)

        LinearAlgebra.$op(n::T1, A::QuantumObject{<:AbstractArray{T2}}) where {T1<:Number,T2} =
            QuantumObject($(op)(n * I, A.data), A.type, A.dims)
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T1}}, n::T2) where {T1,T2<:Number} =
            QuantumObject($(op)(A.data, n * I), A.type, A.dims)
    end
end

function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Ket, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Bra, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},BraQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(vec2mat(A.data * mat2vec(B.data)), Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorBraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorKetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorKetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorKet, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorBraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorBra, A.dims)
end

LinearAlgebra.:(^)(A::QuantumObject{<:AbstractArray{T}}, n::T1) where {T,T1<:Number} =
    QuantumObject(^(A.data, n), A.type, A.dims)
LinearAlgebra.:(/)(A::QuantumObject{<:AbstractArray{T}}, n::T1) where {T,T1<:Number} =
    QuantumObject(/(A.data, n), A.type, A.dims)

@doc raw"""
    dot(A::QuantumObject, B::QuantumObject)

Compute the dot product between two [`QuantumObject`]: ``\langle A | B \rangle``

Note that `A` and `B` should be [`Ket`](@ref) or [`OperatorKet`](@ref)

`A ⋅ B` (where `⋅` can be typed by tab-completing `\\cdot` in the REPL) is a synonym for `dot(A, B)`
"""
function LinearAlgebra.dot(
    A::QuantumObject{<:AbstractArray{T1},OpType},
    B::QuantumObject{<:AbstractArray{T2},OpType},
) where {T1<:Number,T2<:Number,OpType<:Union{KetQuantumObject,OperatorKetQuantumObject}}
    A.dims != B.dims && throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(A.data, B.data)
end

@doc raw"""
    dot(i::QuantumObject, A::QuantumObject j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between three [`QuantumObject`](@ref): ``\langle i | A | j \rangle``

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)
"""
function LinearAlgebra.dot(
    i::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},KetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number}
    ((i.dims != A.dims) || (A.dims != j.dims)) &&
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(i.data, A.data, j.data)
end
function LinearAlgebra.dot(
    i::QuantumObject{<:AbstractArray{T1},OperatorKetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},OperatorKetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number}
    ((i.dims != A.dims) || (A.dims != j.dims)) &&
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(i.data, A.data, j.data)
end

@doc raw"""
    matrix_element(i::QuantumObject, A::QuantumObject j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between three [`QuantumObject`](@ref): ``\langle i | A | j \rangle``

Note that this function is same as `dot(i, A, j)`

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)
"""
matrix_element(
    i::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},KetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number} = dot(i, A, j)
matrix_element(
    i::QuantumObject{<:AbstractArray{T1},OperatorKetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},OperatorKetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number} = dot(i, A, j)

@doc raw"""
    conj(A::QuantumObject)

Return the element-wise complex conjugation of the [`QuantumObject`](@ref).
"""
Base.conj(A::QuantumObject{<:AbstractArray{T}}) where {T} = QuantumObject(conj(A.data), A.type, A.dims)

@doc raw"""
    transpose(A::QuantumObject)

Lazy matrix transpose of the [`QuantumObject`](@ref).
"""
LinearAlgebra.transpose(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(transpose(A.data), A.type, A.dims)

@doc raw"""
    trans(A::QuantumObject)

Lazy matrix transpose of the [`QuantumObject`](@ref).

Note that this function is same as `transpose(A)`
"""
trans(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = transpose(A)

@doc raw"""
    A'
    adjoint(A::QuantumObject)

Lazy adjoint (conjugate transposition) of the [`QuantumObject`](@ref)

Note that `A'` is a synonym for `adjoint(A)`
"""
LinearAlgebra.adjoint(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(adjoint(A.data), A.type, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), Bra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},BraQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), Ket, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},OperatorKetQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), OperatorBra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},OperatorBraQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), OperatorKet, A.dims)

@doc raw"""
    dag(A::QuantumObject)

Lazy adjoint (conjugate transposition) of the [`QuantumObject`](@ref)

Note that this function is same as `adjoint(A)`
"""
dag(A::QuantumObject{<:AbstractArray{T}}) where {T} = adjoint(A)

@doc raw"""
    dagger(A::QuantumObject)

Lazy adjoint (conjugate transposition) of the [`QuantumObject`](@ref)

Note that this function is same as `adjoint(A)`
"""
dagger(A::QuantumObject{<:AbstractArray{T}}) where {T} = adjoint(A)

LinearAlgebra.inv(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(sparse(inv(Matrix(A.data))), A.type, A.dims)

LinearAlgebra.Hermitian(
    A::QuantumObject{<:AbstractArray{T},OpType},
    uplo::Symbol = :U,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dims)

@doc raw"""
    tr(A::QuantumObject})

Returns the trace of `A`.

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

julia> tr(a' * a)
190.0 + 0.0im
```
"""
LinearAlgebra.tr(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    ishermitian(A) ? real(tr(A.data)) : tr(A.data)

@doc raw"""
    svdvals(A::QuantumObject)

Return the singular values of a [`QuantumObject`](@ref) in descending order
"""
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractVector}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:DenseMatrix}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractSparseMatrix}) = svdvals(sparse_to_dense(A.data))

@doc raw"""
    norm(A::QuantumObject, p::Real=2)

If `A` is either [`Ket`](@ref), [`Bra`](@ref), [`OperatorKet`](@ref), or [`OperatorBra`](@ref), returns the standard vector `p`-norm of `A`.
If `A` is either [`Operator`](@ref) or [`SuperOperator`](@ref), returns [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm of `A`.

Note that the default value of `p=2`

# Examples

```
julia> ψ = fock(10, 2)
Quantum Object:   type=Ket   dims=[10]   size=(10,)
10-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> norm(ψ)
1.0
```
"""
LinearAlgebra.norm(
    A::QuantumObject{<:AbstractArray{T},OpType},
    p::Real = 2,
) where {T,OpType<:Union{KetQuantumObject,BraQuantumObject,OperatorKetQuantumObject,OperatorBraQuantumObject}} =
    norm(A.data, p)
function LinearAlgebra.norm(
    A::QuantumObject{<:AbstractArray{T},OpType},
    p::Real = 2,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    p == 2.0 && return norm(A.data, 2)
    return norm(svdvals(A), p)
end
LinearAlgebra.normalize(A::QuantumObject{<:AbstractArray{T}}) where {T} =
    QuantumObject(normalize(A.data), A.type, A.dims)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T}}) where {T} = (normalize!(A.data); A)

LinearAlgebra.triu!(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (triu!(A.data, k); A)
LinearAlgebra.tril!(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (tril!(A.data, k); A)
LinearAlgebra.triu(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(triu(A.data, k), A.type, A.dims)
LinearAlgebra.tril(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(tril(A.data, k), A.type, A.dims)

LinearAlgebra.lmul!(a::Number, B::QuantumObject{<:AbstractArray}) = (lmul!(a, B.data); B)
LinearAlgebra.rmul!(B::QuantumObject{<:AbstractArray}, a::Number) = (rmul!(B.data, a); B)

@inline LinearAlgebra.mul!(y::AbstractVector{Ty}, A::QuantumObject{<:AbstractMatrix{Ta}}, x, α, β) where {Ty,Ta} =
    mul!(y, A.data, x, α, β)

LinearAlgebra.sqrt(A::QuantumObject{<:AbstractArray{T}}) where {T} =
    QuantumObject(sqrt(sparse_to_dense(A.data)), A.type, A.dims)

@doc raw"""
    sqrtm(A::QuantumObject)

Matrix square root of [`Operator`](@ref) type of [`QuantumObject`](@ref)

Note that for other types of [`QuantumObject`](@ref) use `sprt(A)` instead.
"""
sqrtm(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = sqrt(A)

LinearAlgebra.exp(A::QuantumObject{<:AbstractMatrix{T}}) where {T} =
    QuantumObject(dense_to_sparse(exp(A.data)), A.type, A.dims)
LinearAlgebra.exp(A::QuantumObject{<:AbstractSparseMatrix{T}}) where {T} = QuantumObject(_spexp(A.data), A.type, A.dims)

function _spexp(A::SparseMatrixCSC{T,M}; threshold = 1e-14, nonzero_tol = 1e-20) where {T,M}
    m = checksquare(A) # Throws exception if not square

    mat_norm = norm(A, Inf)
    mat_norm == 0 && return eye(m).data
    scaling_factor = nextpow(2, mat_norm) # Native routine, faster
    A = A ./ scaling_factor
    delta = 1

    P = spdiagm(0 => ones(T, m))
    next_term = P
    n = 1

    while delta > threshold
        next_term *= A / n
        if nnz(next_term) / length(next_term) > 0.25
            tidyup!(next_term, nonzero_tol)
        end
        delta = norm(next_term, Inf)
        P += next_term
        n += 1
    end
    for n in 1:log2(scaling_factor)
        P = P * P
        if nnz(P) / length(P) > 0.25
            tidyup!(P, nonzero_tol)
        end
    end
    return P
end

@doc raw"""
    sinm(O::QuantumObject)

Generates the sine of the operator `O`, defined as

``\sin \left( \hat{O} \right) = \frac{e^{i \hat{O}} - e^{-i \hat{O}}}{2 i}``
"""
sinm(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = -0.5im * (exp(1im * O) - exp(-1im * O))

@doc raw"""
    cosm(O::QuantumObject)

Generates the cosine of the operator `O`, defined as

``\cos \left( \hat{O} \right) = \frac{e^{i \hat{O}} + e^{-i \hat{O}}}{2}``
"""
cosm(O::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = 0.5 * (exp(1im * O) + exp(-1im * O))

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

@doc raw"""
    tidyup(A::QuantumObject, tol::Real=1e-14)

Removes those elements of a QuantumObject `A` whose absolute value is less than `tol`.
"""
tidyup(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} =
    QuantumObject(tidyup(A.data, tol), A.type, A.dims)
tidyup(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} = @. T(abs(A) > tol) * A
tidyup(A::AbstractSparseMatrix{T}, tol::T2 = 1e-14) where {T,T2<:Real} = droptol!(copy(A), tol)

@doc raw"""
    tidyup!(A::QuantumObject, tol::Real=1e-14)

Removes those elements of a QuantumObject `A` whose absolute value is less than `tol`.
In-place version of [`tidyup`](#tidyup).
"""
tidyup!(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} = (tidyup!(A.data, tol); A)
tidyup!(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} = @. A = T(abs(A) > tol) * A
tidyup!(A::AbstractSparseMatrix{T}, tol::T2 = 1e-14) where {T,T2<:Real} = droptol!(A, tol)

@doc raw"""
    get_data(A::QuantumObject)

Returns the data of a QuantumObject.
"""
get_data(A::QuantumObject) = A.data

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
