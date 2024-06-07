#=
Arithmetic and Attributes for QuantumObject
    - extend most of the useful functions in LinearAlgebra for QuantumObject
    - export most of the attribute functions in "Python Qobj class"
=#

export trans, dag, dagger, matrix_element, unit
export sqrtm, logm, expm, sinm, cosm
export proj, ptrace, purity
export tidyup, tidyup!
export get_data, get_coherence
export permute

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
trans(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = transpose(A)

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

@doc raw"""
    inv(A::QuantumObject)

Matrix inverse of the [`QuantumObject`](@ref)
"""
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
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractMatrix}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractSparseMatrix}) = svdvals(sparse_to_dense(A.data))

@doc raw"""
    norm(A::QuantumObject, p::Real)

If `A` is either [`Ket`](@ref), [`Bra`](@ref), [`OperatorKet`](@ref), or [`OperatorBra`](@ref), returns the standard vector `p`-norm (default `p=2`) of `A`.
If `A` is either [`Operator`](@ref) or [`SuperOperator`](@ref), returns [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm (default `p=1`) of `A`.

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
    p::Real = 1,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    p == 2.0 && return norm(A.data, 2)
    return norm(svdvals(A), p)
end

@doc raw"""
    normalize(A::QuantumObject, p::Real)

Return normalized [`QuantumObject`](@ref) so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
LinearAlgebra.normalize(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    p::Real = 2,
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} = QuantumObject(A.data / norm(A, p), A.type, A.dims)
LinearAlgebra.normalize(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, p::Real = 1) where {T} =
    QuantumObject(A.data / norm(A, p), A.type, A.dims)

@doc raw"""
    normalize!(A::QuantumObject, p::Real)

Normalize [`QuantumObject`](@ref) in-place so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
LinearAlgebra.normalize!(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    p::Real = 2,
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} = (rmul!(A.data, 1 / norm(A)); A)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, p::Real = 1) where {T} =
    (rmul!(A.data, 1 / norm(A)); A)

@doc raw"""
    unit(A::QuantumObject, p::Real)

Return normalized [`QuantumObject`](@ref) so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

Note that this function is same as `normalize(A, p)`

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
unit(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    p::Real = 2,
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} = normalize(A, p)
unit(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, p::Real = 1) where {T} = normalize(A, p)

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

@doc raw"""
    sqrt(A::QuantumObject)

Square root of [`QuantumObject`](@ref)
"""
LinearAlgebra.sqrt(A::QuantumObject{<:AbstractArray{T}}) where {T} =
    QuantumObject(sqrt(sparse_to_dense(A.data)), A.type, A.dims)

@doc raw"""
    sqrtm(A::QuantumObject)

Matrix square root of [`Operator`](@ref) type of [`QuantumObject`](@ref)

Note that for other types of [`QuantumObject`](@ref) use `sprt(A)` instead.
"""
sqrtm(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = sqrt(A)

@doc raw"""
    log(A::QuantumObject)

Matrix logarithm of [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.log(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(log(sparse_to_dense(A.data)), A.type, A.dims)

@doc raw"""
    logm(A::QuantumObject)

Matrix logarithm of [`QuantumObject`](@ref)

Note that this function is same as `log(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
logm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = log(A)

@doc raw"""
    exp(A::QuantumObject)

Matrix exponential of [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.exp(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(dense_to_sparse(exp(A.data)), A.type, A.dims)
LinearAlgebra.exp(
    A::QuantumObject{<:AbstractSparseMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(_spexp(A.data), A.type, A.dims)

@doc raw"""
    expm(A::QuantumObject)

Matrix exponential of [`QuantumObject`](@ref)

Note that this function is same as `exp(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
expm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = exp(A)

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

Generates the matrix sine of the operator `O`, defined as

``\sin \left( \hat{O} \right) = \frac{e^{i \hat{O}} - e^{-i \hat{O}}}{2 i}``

Note that this function only supports for [`Operator`](@ref)
"""
sinm(
    O::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (exp(1im * O) - exp(-1im * O)) / 2im

@doc raw"""
    cosm(O::QuantumObject)

Generates the matrix cosine of the operator `O`, defined as

``\cos \left( \hat{O} \right) = \frac{e^{i \hat{O}} + e^{-i \hat{O}}}{2}``

Note that this function only supports for [`Operator`](@ref)
"""
cosm(
    O::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (exp(1im * O) + exp(-1im * O)) / 2

@doc raw"""
    diag(A::QuantumObject, k::Int=0)

Return the `k`-th diagonal elements of a matrix-type [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.diag(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
    k::Int = 0,
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = diag(A.data, k)

@doc raw"""
    proj(ψ::QuantumObject)

Return the projector for a [`Ket`](@ref) or [`Bra`](@ref) type of [`QuantumObject`](@ref)
"""
proj(ψ::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} = ψ * ψ'
proj(ψ::QuantumObject{<:AbstractArray{T},BraQuantumObject}) where {T} = ψ' * ψ

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
    purity(ρ::QuantumObject)

Calculate the purity of a [`QuantumObject`](@ref): ``\textrm{Tr}(\rho^2)``

Note that this function only supports for [`Ket`](@ref), [`Bra`](@ref), and [`Operator`](@ref)
"""
purity(ρ::QuantumObject{<:AbstractArray{T},ObjType}) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} =
    sum(abs2, ρ.data)
purity(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = real(tr(ρ.data^2))

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
In-place version of [`tidyup`](@ref).
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

@doc raw"""
    permute(A::QuantumObject, order::Vector{Int})

Permute the tensor structure of a [`QuantumObject`](@ref) `A` according to the specified `order` list

Note that this method currently works for [`Ket`](@ref), [`Bra`](@ref), and [`Operator`](@ref) types of [`QuantumObject`](@ref).

# Examples

If `order = [2, 1, 3]`, the Hilbert space structure will be re-arranged: H₁ ⊗ H₂ ⊗ H₃ → H₂ ⊗ H₁ ⊗ H₃.

```
julia> ψ1 = fock(2, 0)
julia> ψ2 = fock(3, 1)
julia> ψ3 = fock(4, 2)
julia> ψ_123 = tensor(ψ1, ψ2, ψ3)
julia> permute(ψ_123, [2, 1, 3]) ≈ tensor(ψ2, ψ1, ψ3)
true
```
"""
function permute(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    order::AbstractVector{Int},
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject,OperatorQuantumObject}}
    (length(order) != length(A.dims)) &&
        throw(ArgumentError("The order list must have the same length as the number of subsystems (A.dims)"))

    !isperm(order) && throw(ArgumentError("$(order) is not a valid permutation of the subsystems (A.dims)"))

    # obtain the arguments: dims for reshape; perm for PermutedDimsArray
    dims, perm = _dims_and_perm(A.type, A.dims, order, length(order))

    return QuantumObject(reshape(PermutedDimsArray(reshape(A.data, dims...), perm), size(A)), A.type, A.dims[order])
end

function _dims_and_perm(
    ::ObjType,
    dims::AbstractVector{Int},
    order::AbstractVector{Int},
    L::Int,
) where {ObjType<:Union{KetQuantumObject,BraQuantumObject}}
    return reverse(dims), reverse!((L + 1) .- order)
end

function _dims_and_perm(::OperatorQuantumObject, dims::AbstractVector{Int}, order::AbstractVector{Int}, L::Int)
    return reverse!([dims; dims]), reverse!((2 * L + 1) .- [order; order .+ L])
end
