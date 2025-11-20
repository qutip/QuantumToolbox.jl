#=
Arithmetic and Attributes for QuantumObject
    - extend most of the useful functions in LinearAlgebra for QuantumObject
    - export most of the attribute functions in "Python Qobj class"
=#

export proj, ptrace, purity
export tidyup, tidyup!
export get_data, get_coherence

#    Broadcasting
Base.broadcastable(x::QuantumObject) = x.data
for op in (:(+), :(-), :(*), :(/), :(^))
    @eval begin
        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::QuantumObject)
            return QuantumObject(broadcast($op, x.data, y.data), x.type, x.dimensions)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::Number)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dimensions)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::Number, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dimensions)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::AbstractArray)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dimensions)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::AbstractArray, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dimensions)
        end
    end
end

for op in (:(+), :(-), :(*))
    @eval begin
        function Base.$op(A::AbstractQuantumObject, B::AbstractQuantumObject)
            check_dimensions(A, B)
            QType = promote_op_type(A, B)
            return QType($(op)(A.data, B.data), A.type, A.dimensions)
        end
        Base.$op(A::AbstractQuantumObject) = get_typename_wrapper(A)($(op)(A.data), A.type, A.dimensions)

        Base.$op(n::T, A::AbstractQuantumObject) where {T <: Number} =
            get_typename_wrapper(A)($(op)(n * I, A.data), A.type, A.dimensions)
        Base.$op(A::AbstractQuantumObject, n::T) where {T <: Number} =
            get_typename_wrapper(A)($(op)(A.data, n * I), A.type, A.dimensions)
    end
end

function check_mul_dimensions(from::NTuple{NA, AbstractSpace}, to::NTuple{NB, AbstractSpace}) where {NA, NB}
    (from != to) && throw(
        DimensionMismatch(
            "The quantum object with (right) dims = $(dimensions_to_dims(from)) can not multiply a quantum object with (left) dims = $(dimensions_to_dims(to)) on the right-hand side.",
        ),
    )
    return nothing
end

for ADimType in (:ProductDimensions, :GeneralProductDimensions)
    for BDimType in (:ProductDimensions, :GeneralProductDimensions)
        if ADimType == BDimType == :ProductDimensions
            @eval begin
                function Base.:(*)(
                        A::AbstractQuantumObject{Operator, <:$ADimType},
                        B::AbstractQuantumObject{Operator, <:$BDimType},
                    )
                    check_dimensions(A, B)
                    QType = promote_op_type(A, B)
                    return QType(A.data * B.data, Operator(), A.dimensions)
                end
            end
        else
            @eval begin
                function Base.:(*)(
                        A::AbstractQuantumObject{Operator, <:$ADimType},
                        B::AbstractQuantumObject{Operator, <:$BDimType},
                    )
                    check_mul_dimensions(get_dimensions_from(A), get_dimensions_to(B))
                    QType = promote_op_type(A, B)
                    return QType(
                        A.data * B.data,
                        Operator(),
                        GeneralProductDimensions(get_dimensions_to(A), get_dimensions_from(B)),
                    )
                end
            end
        end
    end
end

function Base.:(*)(A::AbstractQuantumObject{Operator}, B::QuantumObject{Ket, <:ProductDimensions})
    check_mul_dimensions(get_dimensions_from(A), get_dimensions_to(B))
    return QuantumObject(A.data * B.data, Ket(), ProductDimensions(get_dimensions_to(A)))
end
function Base.:(*)(A::QuantumObject{Bra, <:ProductDimensions}, B::AbstractQuantumObject{Operator})
    check_mul_dimensions(get_dimensions_from(A), get_dimensions_to(B))
    return QuantumObject(A.data * B.data, Bra(), ProductDimensions(get_dimensions_from(B)))
end
function Base.:(*)(A::QuantumObject{Ket}, B::QuantumObject{Bra})
    check_dimensions(A, B)
    return QuantumObject(A.data * B.data, Operator(), A.dimensions) # to align with QuTiP, don't use kron(A, B) to do it.
end
function Base.:(*)(A::QuantumObject{Bra}, B::QuantumObject{Ket})
    check_dimensions(A, B)
    return A.data * B.data
end
function Base.:(*)(A::AbstractQuantumObject{SuperOperator}, B::QuantumObject{Operator})
    check_dimensions(A, B)
    return QuantumObject(vec2mat(A.data * mat2vec(B.data)), Operator(), A.dimensions)
end
function Base.:(*)(A::QuantumObject{OperatorBra}, B::QuantumObject{OperatorKet})
    check_dimensions(A, B)
    return A.data * B.data
end
function Base.:(*)(A::AbstractQuantumObject{SuperOperator}, B::QuantumObject{OperatorKet})
    check_dimensions(A, B)
    return QuantumObject(A.data * B.data, OperatorKet(), A.dimensions)
end
function Base.:(*)(A::QuantumObject{OperatorBra}, B::AbstractQuantumObject{SuperOperator})
    check_dimensions(A, B)
    return QuantumObject(A.data * B.data, OperatorBra(), A.dimensions)
end

Base.:(^)(A::QuantumObject, n::T) where {T <: Number} = QuantumObject(^(A.data, n), A.type, A.dimensions)
Base.:(/)(A::AbstractQuantumObject, n::T) where {T <: Number} = get_typename_wrapper(A)(A.data / n, A.type, A.dimensions)

@doc raw"""
    A ⋅ B
    dot(A::QuantumObject, B::QuantumObject)

Compute the dot product between two [`QuantumObject`](@ref): ``\langle A | B \rangle``

Note that `A` and `B` should be [`Ket`](@ref) or [`OperatorKet`](@ref)

!!! note
    `A ⋅ B` (where `⋅` can be typed by tab-completing `\cdot` in the REPL) is a synonym of `dot(A, B)`.
"""
function LinearAlgebra.dot(A::QuantumObject{OpType}, B::QuantumObject{OpType}) where {OpType <: Union{Ket, OperatorKet}}
    check_dimensions(A, B)
    return LinearAlgebra.dot(A.data, B.data)
end

@doc raw"""
    dot(i::QuantumObject, A::AbstractQuantumObject j::QuantumObject)
    matrix_element(i::QuantumObject, A::AbstractQuantumObject j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between a [`AbstractQuantumObject`](@ref) and two [`QuantumObject`](@ref) (`i` and `j`), namely ``\langle i | \hat{A} | j \rangle``.

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)

!!! note
    `matrix_element(i, A, j)` is a synonym of `dot(i, A, j)`.
"""
function LinearAlgebra.dot(i::QuantumObject{Ket}, A::AbstractQuantumObject{Operator}, j::QuantumObject{Ket})
    check_dimensions(i, A, j)
    return LinearAlgebra.dot(i.data, A.data, j.data)
end
function LinearAlgebra.dot(
        i::QuantumObject{OperatorKet},
        A::AbstractQuantumObject{SuperOperator},
        j::QuantumObject{OperatorKet},
    )
    check_dimensions(i, A, j)
    return LinearAlgebra.dot(i.data, A.data, j.data)
end

@doc raw"""
    zero(A::AbstractQuantumObject)

Return a similar [`AbstractQuantumObject`](@ref) with `dims` and `type` are same as `A`, but `data` is a zero-array.
"""
Base.zero(A::AbstractQuantumObject) = get_typename_wrapper(A)(zero(A.data), A.type, A.dimensions)

@doc raw"""
    one(A::AbstractQuantumObject)

Return a similar [`AbstractQuantumObject`](@ref) with `dims` and `type` are same as `A`, but `data` is an identity matrix.

Note that `A` must be [`Operator`](@ref) or [`SuperOperator`](@ref).
"""
Base.one(A::AbstractQuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}} =
    get_typename_wrapper(A)(one(A.data), A.type, A.dimensions)

@doc raw"""
    conj(A::AbstractQuantumObject)

Return the element-wise complex conjugation of the [`AbstractQuantumObject`](@ref).
"""
Base.conj(A::AbstractQuantumObject) = get_typename_wrapper(A)(conj(A.data), A.type, A.dimensions)

@doc raw"""
    transpose(A::AbstractQuantumObject)

Lazy matrix transpose of the [`AbstractQuantumObject`](@ref).
"""
Base.transpose(A::AbstractQuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}} =
    get_typename_wrapper(A)(transpose(A.data), A.type, transpose(A.dimensions))

@doc raw"""
    A'
    adjoint(A::AbstractQuantumObject)
    dag(A::AbstractQuantumObject)

Lazy adjoint (conjugate transposition) of the [`AbstractQuantumObject`](@ref)

!!! note
    `A'` and `dag(A)` are synonyms of `adjoint(A)`.
"""
Base.adjoint(A::AbstractQuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}} =
    get_typename_wrapper(A)(adjoint(A.data), A.type, adjoint(A.dimensions))
Base.adjoint(A::QuantumObject{Ket}) = QuantumObject(adjoint(A.data), Bra(), adjoint(A.dimensions))
Base.adjoint(A::QuantumObject{Bra}) = QuantumObject(adjoint(A.data), Ket(), adjoint(A.dimensions))
Base.adjoint(A::QuantumObject{OperatorKet}) = QuantumObject(adjoint(A.data), OperatorBra(), adjoint(A.dimensions))
Base.adjoint(A::QuantumObject{OperatorBra}) = QuantumObject(adjoint(A.data), OperatorKet(), adjoint(A.dimensions))

@doc raw"""
    inv(A::AbstractQuantumObject)

Matrix inverse of the [`AbstractQuantumObject`](@ref). If `A` is a [`QuantumObjectEvolution`](@ref), the inverse is computed at the last computed time.
"""
LinearAlgebra.inv(A::AbstractQuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}} =
    QuantumObject(sparse(inv(Matrix(A.data))), A.type, A.dimensions)

LinearAlgebra.Hermitian(A::QuantumObject{OpType}, uplo::Symbol = :U) where {OpType <: Union{Operator, SuperOperator}} =
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dimensions)

@doc raw"""
    tr(A::QuantumObject)

Returns the trace of [`QuantumObject`](@ref).

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)

# Examples

```jldoctest
julia> a = destroy(20)

Quantum Object:   type=Operator()   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⎡⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⎦

julia> tr(a' * a)
190.0 + 0.0im
```
"""
LinearAlgebra.tr(A::QuantumObject{OpType}) where {OpType <: Union{Operator, SuperOperator}} = tr(A.data)
LinearAlgebra.tr(
    A::QuantumObject{OpType, DimsType, <:Union{<:Hermitian{TF}, Symmetric{TR}}},
) where {OpType <: Operator, DimsType, TF <: Number, TR <: Real} = real(tr(A.data))

@doc raw"""
    svdvals(A::QuantumObject)

Return the singular values of a [`QuantumObject`](@ref) in descending order
"""
LinearAlgebra.svdvals(A::QuantumObject) = svdvals(to_dense(A.data))

@doc raw"""
    norm(A::QuantumObject, p::Real)

Return the standard vector `p`-norm or [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm of a [`QuantumObject`](@ref) depending on the type of `A`:

- If `A` is either [`Ket`](@ref), [`Bra`](@ref), [`OperatorKet`](@ref), or [`OperatorBra`](@ref), returns the standard vector `p`-norm (default `p=2`) of `A`.
- If `A` is either [`Operator`](@ref) or [`SuperOperator`](@ref), returns [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm (default `p=1`) of `A`.

# Examples

```jldoctest
julia> ψ = fock(10, 2)

Quantum Object:   type=Ket()   dims=[10]   size=(10,)
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
LinearAlgebra.norm(A::QuantumObject{OpType}, p::Real = 2) where {OpType <: Union{Ket, Bra, OperatorKet, OperatorBra}} =
    norm(A.data, p)
function LinearAlgebra.norm(A::QuantumObject{OpType}, p::Real = 1) where {OpType <: Union{Operator, SuperOperator}}
    p == 2.0 && return norm(A.data, 2)
    return norm(svdvals(A), p)
end

@doc raw"""
    normalize(A::QuantumObject, p::Real)
    unit(A::QuantumObject, p::Real)

Return normalized [`QuantumObject`](@ref) so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

!!! note
    `unit` is a synonym of `normalize`.

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
LinearAlgebra.normalize(A::QuantumObject{ObjType}, p::Real = 2) where {ObjType <: Union{Ket, Bra}} =
    QuantumObject(A.data / norm(A, p), A.type, A.dimensions)
LinearAlgebra.normalize(A::QuantumObject{Operator}, p::Real = 1) =
    QuantumObject(A.data / norm(A, p), A.type, A.dimensions)

@doc raw"""
    normalize!(A::QuantumObject, p::Real)

Normalize [`QuantumObject`](@ref) in-place so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
LinearAlgebra.normalize!(A::QuantumObject{ObjType}, p::Real = 2) where {ObjType <: Union{Ket, Bra}} =
    (rmul!(A.data, 1 / norm(A, p)); A)
LinearAlgebra.normalize!(A::QuantumObject{Operator}, p::Real = 1) = (rmul!(A.data, 1 / norm(A, p)); A)

LinearAlgebra.triu!(A::QuantumObject{OpType}, k::Integer = 0) where {OpType <: Union{Operator, SuperOperator}} =
    (triu!(A.data, k); A)
LinearAlgebra.tril!(A::QuantumObject{OpType}, k::Integer = 0) where {OpType <: Union{Operator, SuperOperator}} =
    (tril!(A.data, k); A)
LinearAlgebra.triu(A::QuantumObject{OpType}, k::Integer = 0) where {OpType <: Union{Operator, SuperOperator}} =
    QuantumObject(triu(A.data, k), A.type, A.dimensions)
LinearAlgebra.tril(A::QuantumObject{OpType}, k::Integer = 0) where {OpType <: Union{Operator, SuperOperator}} =
    QuantumObject(tril(A.data, k), A.type, A.dimensions)

LinearAlgebra.lmul!(a::Number, B::QuantumObject) = (lmul!(a, B.data); B)
LinearAlgebra.rmul!(B::QuantumObject, a::Number) = (rmul!(B.data, a); B)

@inline LinearAlgebra.mul!(y::AbstractVector{T}, A::QuantumObject, x, α, β) where {T} = mul!(y, A.data, x, α, β)

@doc raw"""
    √(A)
    sqrt(A::QuantumObject)

Matrix square root of [`QuantumObject`](@ref)

!!! note
    `√(A)` (where `√` can be typed by tab-completing `\sqrt` in the REPL) is a synonym of `sqrt(A)`.
"""
Base.sqrt(A::QuantumObject) = QuantumObject(sqrt(to_dense(A.data)), A.type, A.dimensions)

@doc raw"""
    log(A::QuantumObject)

Matrix logarithm of [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
Base.log(A::QuantumObject{ObjType}) where {ObjType <: Union{Operator, SuperOperator}} =
    QuantumObject(log(to_dense(A.data)), A.type, A.dimensions)

@doc raw"""
    exp(A::QuantumObject)

Matrix exponential of [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
Base.exp(A::QuantumObject{ObjType, DimsType, <:AbstractMatrix}) where {ObjType <: Union{Operator, SuperOperator}, DimsType} =
    QuantumObject(to_sparse(exp(A.data)), A.type, A.dimensions)
Base.exp(
    A::QuantumObject{ObjType, DimsType, <:AbstractSparseMatrix},
) where {ObjType <: Union{Operator, SuperOperator}, DimsType} = QuantumObject(_spexp(A.data), A.type, A.dimensions)

function _spexp(A::SparseMatrixCSC{T, M}; threshold = 1.0e-14, nonzero_tol = 1.0e-20) where {T <: Number, M <: Int}
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
Base.sin(A::QuantumObject{ObjType}) where {ObjType <: Union{Operator, SuperOperator}} =
    (exp(1im * A) - exp(-1im * A)) / 2im

@doc raw"""
    cos(A::QuantumObject)

Matrix cosine of [`QuantumObject`](@ref), defined as

``\cos \left( \hat{A} \right) = \frac{e^{i \hat{A}} + e^{-i \hat{A}}}{2}``

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
Base.cos(A::QuantumObject{ObjType}) where {ObjType <: Union{Operator, SuperOperator}} = (exp(1im * A) + exp(-1im * A)) / 2

@doc raw"""
    diag(A::QuantumObject, k::Int=0)

Return the `k`-th diagonal elements of a matrix-type [`QuantumObject`](@ref)

Note that this function only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
LinearAlgebra.diag(A::QuantumObject{ObjType}, k::Int = 0) where {ObjType <: Union{Operator, SuperOperator}} =
    diag(A.data, k)

@doc raw"""
    proj(ψ::QuantumObject)

Return the projector for a [`Ket`](@ref) or [`Bra`](@ref) type of [`QuantumObject`](@ref)
"""
proj(ψ::QuantumObject{Ket}) = ψ * ψ'
proj(ψ::QuantumObject{Bra}) = ψ' * ψ

@doc raw"""
    ptrace(QO::QuantumObject, sel)

[Partial trace](https://en.wikipedia.org/wiki/Partial_trace) of a quantum state `QO` leaving only the dimensions with the indices present in the `sel` vector.

Note that this function will always return [`Operator`](@ref). No matter the input [`QuantumObject`](@ref) is a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).

# Examples

Two qubits in the state ``\ket{\psi} = \ket{e,g}``:
```jldoctest
julia> ψ = kron(fock(2,0), fock(2,1))

Quantum Object:   type=Ket()   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> ptrace(ψ, 2)

Quantum Object:   type=Operator()   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im
```

or in an entangled state ``\ket{\psi} = \frac{1}{\sqrt{2}} \left( \ket{e,e} + \ket{g,g} \right)``:
```jldoctest
julia> ψ = 1 / √2 * (kron(fock(2,0), fock(2,0)) + kron(fock(2,1), fock(2,1)))

Quantum Object:   type=Ket()   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im

julia> ptrace(ψ, 1)

Quantum Object:   type=Operator()   dims=[2]   size=(2, 2)   ishermitian=true
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im
 0.0+0.0im  0.5+0.0im
```
"""
function ptrace(QO::QuantumObject{Ket}, sel::Union{AbstractVector{Int}, Tuple})
    any(s -> s isa EnrSpace, QO.dimensions.to) && throw(ArgumentError("ptrace does not support EnrSpace"))

    _non_static_array_warning("sel", sel)

    if length(sel) == 0 # return full trace for empty sel
        return tr(ket2dm(QO))
    else
        n_d = length(QO.dimensions)

        (any(>(n_d), sel) || any(<(1), sel)) && throw(
            ArgumentError("Invalid indices in `sel`: $(sel), the given QuantumObject only have $(n_d) sub-systems"),
        )
        allunique(sel) || throw(ArgumentError("Duplicate selection indices in `sel`: $(sel)"))
        (n_d == 1) && return ket2dm(QO) # ptrace should always return Operator
    end

    _sort_sel = sort(SVector{length(sel), Int}(sel))
    ρtr, dkeep = _ptrace_ket(QO.data, QO.dims, _sort_sel)
    return QuantumObject(ρtr, type = Operator(), dims = ProductDimensions(dkeep))
end

ptrace(QO::QuantumObject{Bra}, sel::Union{AbstractVector{Int}, Tuple}) = ptrace(QO', sel)

function ptrace(QO::QuantumObject{Operator}, sel::Union{AbstractVector{Int}, Tuple})
    any(s -> s isa EnrSpace, QO.dimensions.to) && throw(ArgumentError("ptrace does not support EnrSpace"))

    # TODO: support for special cases when some of the subsystems have same `to` and `from` space
    isa(QO.dimensions, GeneralProductDimensions) &&
        (get_dimensions_to(QO) != get_dimensions_from(QO)) &&
        throw(ArgumentError("Invalid partial trace for dims = $(_get_dims_string(QO.dimensions))"))

    _non_static_array_warning("sel", sel)

    if length(sel) == 0 # return full trace for empty sel
        return tr(QO)
    else
        n_d = length(QO.dimensions)

        (any(>(n_d), sel) || any(<(1), sel)) && throw(
            ArgumentError("Invalid indices in `sel`: $(sel), the given QuantumObject only have $(n_d) sub-systems"),
        )
        allunique(sel) || throw(ArgumentError("Duplicate selection indices in `sel`: $(sel)"))
        (n_d == 1) && return QO
    end

    dims = dimensions_to_dims(get_dimensions_to(QO))
    _sort_sel = sort(SVector{length(sel), Int}(sel))
    ρtr, dkeep = _ptrace_oper(QO.data, dims, _sort_sel)
    return QuantumObject(ρtr, type = Operator(), dims = ProductDimensions(dkeep))
end
ptrace(QO::QuantumObject, sel::Int) = ptrace(QO, SVector(sel))

function _ptrace_ket(QO::AbstractArray, dims::Union{SVector, MVector}, sel)
    n_d = length(dims)

    n_d == 1 && return QO, dims

    qtrace = filter(i -> i ∉ sel, 1:n_d)
    dkeep = dims[sel]
    dtrace = dims[qtrace]
    n_t = length(dtrace)

    # Concatenate qtrace and sel without losing the length information
    # Tuple(qtrace..., sel...)
    qtrace_sel = ntuple(Val(n_d)) do i
        if i <= n_t
            @inbounds qtrace[i]
        else
            @inbounds sel[i - n_t]
        end
    end

    vmat = reshape(QO, reverse(dims)...)
    topermute = reverse(n_d + 1 .- qtrace_sel)
    vmat = permutedims(vmat, topermute) # TODO: use PermutedDimsArray when Julia v1.11.0 is released
    vmat = reshape(vmat, prod(dkeep), prod(dtrace))

    return vmat * vmat', dkeep
end

function _ptrace_oper(QO::AbstractArray, dims::Union{SVector, MVector}, sel)
    n_d = length(dims)

    n_d == 1 && return QO, dims

    qtrace = filter(i -> i ∉ sel, 1:n_d)
    dkeep = dims[sel]
    dtrace = dims[qtrace]
    n_k = length(dkeep)
    n_t = length(dtrace)
    _2_n_t = 2 * n_t

    # Concatenate qtrace and sel without losing the length information
    # Tuple(qtrace..., sel...)
    qtrace_sel = ntuple(Val(2 * n_d)) do i
        if i <= n_t
            @inbounds qtrace[i]
        elseif i <= _2_n_t
            @inbounds qtrace[i - n_t] + n_d
        elseif i <= _2_n_t + n_k
            @inbounds sel[i - _2_n_t]
        else
            @inbounds sel[i - _2_n_t - n_k] + n_d
        end
    end

    ρmat = reshape(QO, reverse(vcat(dims, dims))...)
    topermute = reverse(2 * n_d + 1 .- qtrace_sel)
    ρmat = permutedims(ρmat, topermute) # TODO: use PermutedDimsArray when Julia v1.11.0 is released
    ρmat = reshape(ρmat, prod(dkeep), prod(dkeep), prod(dtrace), prod(dtrace))
    res = _map_trace(ρmat)

    return res, dkeep
end

_map_trace(A::AbstractArray{T, 4}) where {T} = map(tr, eachslice(A, dims = (1, 2)))

@doc raw"""
    purity(ρ::QuantumObject)

Calculate the purity of a [`QuantumObject`](@ref): ``\textrm{Tr}(\rho^2)``

Note that this function only supports for [`Ket`](@ref), [`Bra`](@ref), and [`Operator`](@ref)
"""
purity(ρ::QuantumObject{ObjType}) where {ObjType <: Union{Ket, Bra}} = sum(abs2, ρ.data)
purity(ρ::QuantumObject{Operator}) = real(tr(ρ.data^2))

@doc raw"""
    tidyup(A::QuantumObject, tol::Real=settings.tidyup_tol)

Given a [`QuantumObject`](@ref) `A`, check the real and imaginary parts of each element separately. Remove the real or imaginary value if its absolute value is less than `tol`.
"""
tidyup(A::QuantumObject, tol::T = settings.tidyup_tol) where {T <: Real} =
    QuantumObject(tidyup(A.data, tol), A.type, A.dimensions)
tidyup(A::AbstractArray, tol::T2 = settings.tidyup_tol) where {T2 <: Real} = tidyup!(copy(A), tol)

@doc raw"""
    tidyup!(A::QuantumObject, tol::Real=settings.tidyup_tol)

Given a [`QuantumObject`](@ref) `A`, check the real and imaginary parts of each element separately. Remove the real or imaginary value if its absolute value is less than `tol`.

Note that this function is an in-place version of [`tidyup`](@ref).
"""
tidyup!(A::QuantumObject, tol::T = settings.tidyup_tol) where {T <: Real} = (tidyup!(A.data, tol); A)
function tidyup!(A::AbstractSparseArray, tol::T2 = settings.tidyup_tol) where {T2 <: Real}
    tidyup!(nonzeros(A), tol) # tidyup A.nzval in-place (also support for CUDA sparse arrays)
    return dropzeros!(A)
end
tidyup!(A::AbstractArray{T}, tol::T2 = settings.tidyup_tol) where {T <: Real, T2 <: Real} = @. A = T(abs(A) > tol) * A
tidyup!(A::AbstractArray{T}, tol::T2 = settings.tidyup_tol) where {T, T2 <: Real} =
    @. A = T(abs(real(A)) > tol) * real(A) + 1im * T(abs(imag(A)) > tol) * imag(A)

@doc raw"""
    get_data(A::AbstractQuantumObject)

Returns the data of a [`AbstractQuantumObject`](@ref).
"""
get_data(A::AbstractQuantumObject) = getfield(A, :data)

@doc raw"""
    get_coherence(ψ::QuantumObject)

Get the coherence value ``\alpha`` by measuring the expectation value of the destruction operator ``\hat{a}`` on a state ``\ket{\psi}`` or a density matrix ``\hat{\rho}``.

It returns both ``\alpha`` and the corresponding state with the coherence removed: ``\ket{\delta_\alpha} = \exp ( \alpha^* \hat{a} - \alpha \hat{a}^\dagger ) \ket{\psi}`` for a pure state, and ``\hat{\rho_\alpha} = \exp ( \alpha^* \hat{a} - \alpha \hat{a}^\dagger ) \hat{\rho} \exp ( -\bar{\alpha} \hat{a} + \alpha \hat{a}^\dagger )`` for a density matrix. These states correspond to the quantum fluctuations around the coherent state ``\ket{\alpha}`` or ``|\alpha\rangle\langle\alpha|``.
"""
function get_coherence(ψ::QuantumObject{Ket})
    a = destroy(hilbert_dimensions_to_size(ψ.dimensions)[1])
    α = expect(a, ψ)
    D = exp(α * a' - conj(α) * a)

    return α, D' * ψ
end

function get_coherence(ρ::QuantumObject{Operator})
    a = destroy(hilbert_dimensions_to_size(ρ.dimensions)[1])
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

```jldoctest
julia> ψ1 = fock(2, 0);

julia> ψ2 = fock(3, 1);

julia> ψ3 = fock(4, 2);

julia> ψ_123 = tensor(ψ1, ψ2, ψ3);

julia> permute(ψ_123, (2, 1, 3)) ≈ tensor(ψ2, ψ1, ψ3)
true
```

!!! warning "Beware of type-stability!"
    It is highly recommended to use `permute(A, order)` with `order` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function SparseArrays.permute(
        A::QuantumObject{ObjType},
        order::Union{AbstractVector{Int}, Tuple},
    ) where {ObjType <: Union{Ket, Bra, Operator}}
    any(s -> s isa EnrSpace, A.dimensions.to) && throw(ArgumentError("permute does not support EnrSpace"))

    (length(order) != length(A.dimensions)) &&
        throw(ArgumentError("The order list must have the same length as the number of subsystems (A.dims)"))

    !isperm(order) && throw(ArgumentError("$(order) is not a valid permutation of the subsystems (A.dims)"))

    _non_static_array_warning("order", order)

    order_svector = SVector{length(order), Int}(order) # convert it to SVector for performance

    # obtain the arguments: dims for reshape; perm for PermutedDimsArray
    dims, perm = _dims_and_perm(A.type, A.dims, order_svector, length(order_svector))

    order_dimensions = _order_dimensions(A.dimensions, order_svector)

    return QuantumObject(reshape(permutedims(reshape(A.data, dims...), Tuple(perm)), size(A)), A.type, order_dimensions)
end

_dims_and_perm(::ObjType, dims::SVector{N, Int}, order::AbstractVector{Int}, L::Int) where {ObjType <: Union{Ket, Bra}, N} =
    reverse(dims), reverse((L + 1) .- order)

# if dims originates from ProductDimensions
_dims_and_perm(::Operator, dims::SVector{N, Int}, order::AbstractVector{Int}, L::Int) where {N} =
    reverse(vcat(dims, dims)), reverse((2 * L + 1) .- vcat(order, order .+ L))

# if dims originates from GeneralProductDimensions
_dims_and_perm(::Operator, dims::SVector{2, SVector{N, Int}}, order::AbstractVector{Int}, L::Int) where {N} =
    reverse(vcat(dims[2], dims[1])), reverse((2 * L + 1) .- vcat(order, order .+ L))

_order_dimensions(dimensions::ProductDimensions, order::AbstractVector{Int}) = ProductDimensions(dimensions.to[order])
_order_dimensions(dimensions::GeneralProductDimensions, order::AbstractVector{Int}) =
    GeneralProductDimensions(dimensions.to[order], dimensions.from[order])
