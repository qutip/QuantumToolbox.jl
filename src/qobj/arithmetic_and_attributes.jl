#=
Arithmetic and Attributes for QuantumObject
    - extend most of the useful functions in LinearAlgebra for QuantumObject
    - export most of the attribute functions in "Python Qobj class"
=#

export proj, ptrace, purity, permute
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
        function LinearAlgebra.$op(A::AbstractQuantumObject, B::AbstractQuantumObject)
            A.dims != B.dims &&
                throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
            if A isa QuantumObjectEvolution || B isa QuantumObjectEvolution
                return QuantumObjectEvolution($(op)(A.data, B.data), A.type, A.dims)
            end
            return QuantumObject($(op)(A.data, B.data), A.type, A.dims)
        end
        LinearAlgebra.$op(A::AbstractQuantumObject) = get_typename_wrapper(A)($(op)(A.data), A.type, A.dims)

        LinearAlgebra.$op(n::T, A::AbstractQuantumObject) where {T<:Number} =
            get_typename_wrapper(A)($(op)(n * I, A.data), A.type, A.dims)
        LinearAlgebra.$op(A::AbstractQuantumObject, n::T) where {T<:Number} =
            get_typename_wrapper(A)($(op)(A.data, n * I), A.type, A.dims)
    end
end

function LinearAlgebra.:(*)(
    A::AbstractQuantumObject{DT1,OperatorQuantumObject},
    B::QuantumObject{DT2,KetQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Ket, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{DT1,BraQuantumObject},
    B::AbstractQuantumObject{DT2,OperatorQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Bra, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{DT1,KetQuantumObject},
    B::QuantumObject{DT2,BraQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{DT1,BraQuantumObject},
    B::QuantumObject{DT2,KetQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::AbstractQuantumObject{DT1,SuperOperatorQuantumObject},
    B::QuantumObject{DT2,OperatorQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(vec2mat(A.data * mat2vec(B.data)), Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{DT1,OperatorBraQuantumObject},
    B::QuantumObject{DT2,OperatorKetQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::AbstractQuantumObject{DT1,SuperOperatorQuantumObject},
    B::QuantumObject{DT2,OperatorKetQuantumObject},
) where {DT1,DT2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorKet, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorBraQuantumObject},
    B::AbstractQuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(DimensionMismatch("The two quantum object don't have the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorBra, A.dims)
end

LinearAlgebra.:(^)(A::QuantumObject{DT}, n::T) where {DT,T<:Number} = QuantumObject(^(A.data, n), A.type, A.dims)
LinearAlgebra.:(/)(A::QuantumObject{DT}, n::T) where {DT,T<:Number} = QuantumObject(/(A.data, n), A.type, A.dims)

@doc raw"""
    dot(A::QuantumObject, B::QuantumObject)

Compute the dot product between two [`QuantumObject`](@ref): ``\langle A | B \rangle``

Note that `A` and `B` should be [`Ket`](@ref) or [`OperatorKet`](@ref)

`A ⋅ B` (where `⋅` can be typed by tab-completing `\cdot` in the REPL) is a synonym for `dot(A, B)`
"""
function LinearAlgebra.dot(
    A::QuantumObject{DT1,OpType},
    B::QuantumObject{DT2,OpType},
) where {DT1,DT2,OpType<:Union{KetQuantumObject,OperatorKetQuantumObject}}
    A.dims != B.dims && throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(A.data, B.data)
end

@doc raw"""
    dot(i::QuantumObject, A::AbstractQuantumObject j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between three [`AbstractQuantumObject`](@ref): ``\langle i | \hat{A} | j \rangle``

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)
"""
function LinearAlgebra.dot(
    i::QuantumObject{DT1,KetQuantumObject},
    A::AbstractQuantumObject{DT2,OperatorQuantumObject},
    j::QuantumObject{DT3,KetQuantumObject},
) where {DT1,DT2,DT3}
    ((i.dims != A.dims) || (A.dims != j.dims)) &&
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(i.data, A.data, j.data)
end
function LinearAlgebra.dot(
    i::QuantumObject{DT1,OperatorKetQuantumObject},
    A::AbstractQuantumObject{DT2,SuperOperatorQuantumObject},
    j::QuantumObject{DT3,OperatorKetQuantumObject},
) where {DT1,DT2,DT3}
    ((i.dims != A.dims) || (A.dims != j.dims)) &&
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(i.data, A.data, j.data)
end

@doc raw"""
    conj(A::AbstractQuantumObject)

Return the element-wise complex conjugation of the [`AbstractQuantumObject`](@ref).
"""
Base.conj(A::AbstractQuantumObject) = get_typename_wrapper(A)(conj(A.data), A.type, A.dims)

@doc raw"""
    transpose(A::AbstractQuantumObject)

Lazy matrix transpose of the [`AbstractQuantumObject`](@ref).
"""
LinearAlgebra.transpose(
    A::AbstractQuantumObject{DT,OpType},
) where {DT,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    get_typename_wrapper(A)(transpose(A.data), A.type, A.dims)

@doc raw"""
    A'
    adjoint(A::AbstractQuantumObject)

Lazy adjoint (conjugate transposition) of the [`AbstractQuantumObject`](@ref)

Note that `A'` is a synonym for `adjoint(A)`
"""
LinearAlgebra.adjoint(
    A::AbstractQuantumObject{DT,OpType},
) where {DT,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    get_typename_wrapper(A)(adjoint(A.data), A.type, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{DT,KetQuantumObject}) where {DT} = QuantumObject(adjoint(A.data), Bra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{DT,BraQuantumObject}) where {DT} = QuantumObject(adjoint(A.data), Ket, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{DT,OperatorKetQuantumObject}) where {DT} =
    QuantumObject(adjoint(A.data), OperatorBra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{DT,OperatorBraQuantumObject}) where {DT} =
    QuantumObject(adjoint(A.data), OperatorKet, A.dims)

@doc raw"""
    inv(A::AbstractQuantumObject)

Matrix inverse of the [`AbstractQuantumObject`](@ref). If `A` is a [`QuantumObjectEvolution`](@ref), the inverse is computed at the last computed time.
"""
LinearAlgebra.inv(
    A::AbstractQuantumObject{DT,OpType},
) where {DT,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(sparse(inv(Matrix(A.data))), A.type, A.dims)

LinearAlgebra.Hermitian(
    A::QuantumObject{DT,OpType},
    uplo::Symbol = :U,
) where {DT,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dims)

@doc raw"""
    tr(A::QuantumObject)

Returns the trace of [`QuantumObject`](@ref).

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)

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
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = tr(A.data)
LinearAlgebra.tr(
    A::QuantumObject{<:Union{<:Hermitian{TF},Symmetric{TR}},OpType},
) where {TF<:BlasFloat,TR<:Real,OpType<:OperatorQuantumObject} = real(tr(A.data))

@doc raw"""
    svdvals(A::QuantumObject)

Return the singular values of a [`QuantumObject`](@ref) in descending order
"""
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractVector}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractMatrix}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractSparseMatrix}) = svdvals(sparse_to_dense(A.data))

@doc raw"""
    norm(A::QuantumObject, p::Real)

Return the standard vector `p`-norm or [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm of a [`QuantumObject`](@ref) depending on the type of `A`:

- If `A` is either [`Ket`](@ref), [`Bra`](@ref), [`OperatorKet`](@ref), or [`OperatorBra`](@ref), returns the standard vector `p`-norm (default `p=2`) of `A`.
- If `A` is either [`Operator`](@ref) or [`SuperOperator`](@ref), returns [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm (default `p=1`) of `A`.

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
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} = (rmul!(A.data, 1 / norm(A, p)); A)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, p::Real = 1) where {T} =
    (rmul!(A.data, 1 / norm(A, p)); A)

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
    √(A)
    sqrt(A::QuantumObject)

Matrix square root of [`QuantumObject`](@ref)

Note that `√(A)` is a synonym for `sqrt(A)`
"""
LinearAlgebra.sqrt(A::QuantumObject{<:AbstractArray{T}}) where {T} =
    QuantumObject(sqrt(sparse_to_dense(A.data)), A.type, A.dims)

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

function _spexp(A::SparseMatrixCSC{T,M}; threshold = 1e-14, nonzero_tol = 1e-20) where {T,M}
    m = checksquare(A) # Throws exception if not square

    mat_norm = norm(A, Inf)
    mat_norm == 0 && return sparse(T(1) * I, m, m)
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
    sin(A::QuantumObject)

Matrix sine of [`QuantumObject`](@ref), defined as

``\sin \left( \hat{A} \right) = \frac{e^{i \hat{A}} - e^{-i \hat{A}}}{2 i}``

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.sin(
    A::QuantumObject{DT,ObjType},
) where {DT,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (exp(1im * A) - exp(-1im * A)) / 2im

@doc raw"""
    cos(A::QuantumObject)

Matrix cosine of [`QuantumObject`](@ref), defined as

``\cos \left( \hat{A} \right) = \frac{e^{i \hat{A}} + e^{-i \hat{A}}}{2}``

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.cos(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (exp(1im * A) + exp(-1im * A)) / 2

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
    ptrace(QO::QuantumObject, sel)

[Partial trace](https://en.wikipedia.org/wiki/Partial_trace) of a quantum state `QO` leaving only the dimensions with the indices present in the `sel` vector.

Note that this function will always return [`Operator`](@ref). No matter the input [`QuantumObject`](@ref) is a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).

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

julia> ptrace(ψ, 2)
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

julia> ptrace(ψ, 1)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im
```
"""
function ptrace(QO::QuantumObject{<:AbstractArray,KetQuantumObject}, sel::Union{AbstractVector{Int},Tuple})
    _non_static_array_warning("sel", sel)

    ns = length(sel)
    if ns == 0 # return full trace for empty sel
        return tr(ket2dm(QO))
    else
        nd = length(QO.dims)

        (any(>(nd), sel) || any(<(1), sel)) && throw(
            ArgumentError("Invalid indices in `sel`: $(sel), the given QuantumObject only have $(nd) sub-systems"),
        )
        (ns != length(unique(sel))) && throw(ArgumentError("Duplicate selection indices in `sel`: $(sel)"))
        (nd == 1) && return ket2dm(QO) # ptrace should always return Operator
    end

    _sort_sel = sort(SVector{length(sel),Int}(sel))
    ρtr, dkeep = _ptrace_ket(QO.data, QO.dims, _sort_sel)
    return QuantumObject(ρtr, type = Operator, dims = dkeep)
end

ptrace(QO::QuantumObject{<:AbstractArray,BraQuantumObject}, sel::Union{AbstractVector{Int},Tuple}) = ptrace(QO', sel)

function ptrace(QO::QuantumObject{<:AbstractArray,OperatorQuantumObject}, sel::Union{AbstractVector{Int},Tuple})
    _non_static_array_warning("sel", sel)

    ns = length(sel)
    if ns == 0 # return full trace for empty sel
        return tr(QO)
    else
        nd = length(QO.dims)

        (any(>(nd), sel) || any(<(1), sel)) && throw(
            ArgumentError("Invalid indices in `sel`: $(sel), the given QuantumObject only have $(nd) sub-systems"),
        )
        (ns != length(unique(sel))) && throw(ArgumentError("Duplicate selection indices in `sel`: $(sel)"))
        (nd == 1) && return QO
    end

    _sort_sel = sort(SVector{length(sel),Int}(sel))
    ρtr, dkeep = _ptrace_oper(QO.data, QO.dims, _sort_sel)
    return QuantumObject(ρtr, type = Operator, dims = dkeep)
end
ptrace(QO::QuantumObject, sel::Int) = ptrace(QO, SVector(sel))

function _ptrace_ket(QO::AbstractArray, dims::Union{SVector,MVector}, sel)
    nd = length(dims)

    nd == 1 && return QO, dims

    qtrace = filter(i -> i ∉ sel, 1:nd)
    dkeep = dims[sel]
    dtrace = dims[qtrace]
    nt = length(dtrace)

    # Concatenate qtrace and sel without losing the length information
    # Tuple(qtrace..., sel...)
    qtrace_sel = ntuple(Val(nd)) do i
        if i <= nt
            @inbounds qtrace[i]
        else
            @inbounds sel[i-nt]
        end
    end

    vmat = reshape(QO, reverse(dims)...)
    topermute = reverse(nd + 1 .- qtrace_sel)
    vmat = permutedims(vmat, topermute) # TODO: use PermutedDimsArray when Julia v1.11.0 is released
    vmat = reshape(vmat, prod(dkeep), prod(dtrace))

    return vmat * vmat', dkeep
end

function _ptrace_oper(QO::AbstractArray, dims::Union{SVector,MVector}, sel)
    nd = length(dims)

    nd == 1 && return QO, dims

    qtrace = filter(i -> i ∉ sel, 1:nd)
    dkeep = dims[sel]
    dtrace = dims[qtrace]
    nk = length(dkeep)
    nt = length(dtrace)
    _2_nt = 2 * nt

    # Concatenate qtrace and sel without losing the length information
    # Tuple(qtrace..., sel...)
    qtrace_sel = ntuple(Val(2 * nd)) do i
        if i <= nt
            @inbounds qtrace[i]
        elseif i <= _2_nt
            @inbounds qtrace[i-nt] + nd
        elseif i <= _2_nt + nk
            @inbounds sel[i-_2_nt]
        else
            @inbounds sel[i-_2_nt-nk] + nd
        end
    end

    ρmat = reshape(QO, reverse(vcat(dims, dims))...)
    topermute = reverse(2 * nd + 1 .- qtrace_sel)
    ρmat = permutedims(ρmat, topermute) # TODO: use PermutedDimsArray when Julia v1.11.0 is released
    ρmat = reshape(ρmat, prod(dkeep), prod(dkeep), prod(dtrace), prod(dtrace))
    res = map(tr, eachslice(ρmat, dims = (1, 2)))

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

Given a [`QuantumObject`](@ref) `A`, check the real and imaginary parts of each element separately. Remove the real or imaginary value if its absolute value is less than `tol`.
"""
tidyup(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} =
    QuantumObject(tidyup(A.data, tol), A.type, A.dims)
tidyup(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} = tidyup!(copy(A), tol)

@doc raw"""
    tidyup!(A::QuantumObject, tol::Real=1e-14)

Given a [`QuantumObject`](@ref) `A`, check the real and imaginary parts of each element separately. Remove the real or imaginary value if its absolute value is less than `tol`.

Note that this function is an in-place version of [`tidyup`](@ref).
"""
tidyup!(A::QuantumObject{<:AbstractArray{T}}, tol::T2 = 1e-14) where {T,T2<:Real} = (tidyup!(A.data, tol); A)
function tidyup!(A::AbstractSparseArray{T}, tol::T2 = 1e-14) where {T,T2<:Real}
    tidyup!(nonzeros(A), tol) # tidyup A.nzval in-place (also support for CUDA sparse arrays)
    return dropzeros!(A)
end
tidyup!(A::AbstractArray{T}, tol::T2 = 1e-14) where {T<:Real,T2<:Real} = @. A = T(abs(A) > tol) * A
tidyup!(A::AbstractArray{T}, tol::T2 = 1e-14) where {T,T2<:Real} =
    @. A = T(abs(real(A)) > tol) * real(A) + 1im * T(abs(imag(A)) > tol) * imag(A)

@doc raw"""
    get_data(A::AbstractQuantumObject)

Returns the data of a [`AbstractQuantumObject`](@ref).
"""
get_data(A::AbstractQuantumObject) = A.data

@doc raw"""
    get_coherence(ψ::QuantumObject)

Get the coherence value ``\alpha`` by measuring the expectation value of the destruction operator ``\hat{a}`` on a state ``\ket{\psi}`` or a density matrix ``\hat{\rho}``.

It returns both ``\alpha`` and the corresponding state with the coherence removed: ``\ket{\delta_\alpha} = \exp ( \alpha^* \hat{a} - \alpha \hat{a}^\dagger ) \ket{\psi}`` for a pure state, and ``\hat{\rho_\alpha} = \exp ( \alpha^* \hat{a} - \alpha \hat{a}^\dagger ) \hat{\rho} \exp ( -\bar{\alpha} \hat{a} + \alpha \hat{a}^\dagger )`` for a density matrix. These states correspond to the quantum fluctuations around the coherent state ``\ket{\alpha}`` or ``|\alpha\rangle\langle\alpha|``.
"""
function get_coherence(ψ::QuantumObject{<:AbstractArray,KetQuantumObject})
    a = destroy(prod(ψ.dims))
    α = expect(a, ψ)
    D = exp(α * a' - conj(α) * a)

    return α, D' * ψ
end

function get_coherence(ρ::QuantumObject{<:AbstractArray,OperatorQuantumObject})
    a = destroy(prod(ρ.dims))
    α = expect(a, ρ)
    D = exp(α * a' - conj(α) * a)

    return α, D' * ρ * D
end

@doc raw"""
    permute(A::QuantumObject, order::Union{AbstractVector{Int},Tuple})

Permute the tensor structure of a [`QuantumObject`](@ref) `A` according to the specified `order` list

Note that this method currently works for [`Ket`](@ref), [`Bra`](@ref), and [`Operator`](@ref) types of [`QuantumObject`](@ref).

# Examples

If `order = [2, 1, 3]`, the Hilbert space structure will be re-arranged: ``\mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \mathcal{H}_3 \rightarrow \mathcal{H}_2 \otimes \mathcal{H}_1 \otimes \mathcal{H}_3``.

```
julia> ψ1 = fock(2, 0)
julia> ψ2 = fock(3, 1)
julia> ψ3 = fock(4, 2)
julia> ψ_123 = tensor(ψ1, ψ2, ψ3)
julia> permute(ψ_123, [2, 1, 3]) ≈ tensor(ψ2, ψ1, ψ3)
true
```

!!! warning "Beware of type-stability!"
    It is highly recommended to use `permute(A, order)` with `order` as `Tuple` or `SVector` to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function permute(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    order::Union{AbstractVector{Int},Tuple},
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject,OperatorQuantumObject}}
    (length(order) != length(A.dims)) &&
        throw(ArgumentError("The order list must have the same length as the number of subsystems (A.dims)"))

    !isperm(order) && throw(ArgumentError("$(order) is not a valid permutation of the subsystems (A.dims)"))

    _non_static_array_warning("order", order)

    order_svector = SVector{length(order),Int}(order) # convert it to SVector for performance

    # obtain the arguments: dims for reshape; perm for PermutedDimsArray
    dims, perm = _dims_and_perm(A.type, A.dims, order_svector, length(order_svector))

    return QuantumObject(
        reshape(permutedims(reshape(A.data, dims...), Tuple(perm)), size(A)),
        A.type,
        A.dims[order_svector],
    )
end

function _dims_and_perm(
    ::ObjType,
    dims::SVector{N,Int},
    order::AbstractVector{Int},
    L::Int,
) where {ObjType<:Union{KetQuantumObject,BraQuantumObject},N}
    return reverse(dims), reverse((L + 1) .- order)
end

function _dims_and_perm(::OperatorQuantumObject, dims::SVector{N,Int}, order::AbstractVector{Int}, L::Int) where {N}
    return reverse(vcat(dims, dims)), reverse((2 * L + 1) .- vcat(order, order .+ L))
end
