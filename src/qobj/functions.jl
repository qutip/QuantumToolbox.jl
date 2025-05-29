#=
Functions which manipulates QuantumObject
=#

export ket2dm
export expect, variance
export to_dense, to_sparse
export vec2mat, mat2vec

@doc raw"""
    ket2dm(ψ::QuantumObject)

Transform the ket state ``\ket{\psi}`` into a pure density matrix ``\hat{\rho} = |\psi\rangle\langle\psi|``.
"""
ket2dm(ψ::QuantumObject{Ket}) = ψ * ψ'

ket2dm(ρ::QuantumObject{Operator}) = ρ

@doc raw"""
    expect(O::Union{AbstractQuantumObject,Vector{AbstractQuantumObject}}, ψ::Union{QuantumObject,Vector{QuantumObject}})

Expectation value of the [`Operator`](@ref) `O` with the state `ψ`. The state can be a [`Ket`](@ref), [`Bra`](@ref) or [`Operator`](@ref).

If `ψ` is a [`Ket`](@ref) or [`Bra`](@ref), the function calculates ``\langle\psi|\hat{O}|\psi\rangle``.

If `ψ` is a density matrix ([`Operator`](@ref)), the function calculates ``\textrm{Tr} \left[ \hat{O} \hat{\psi} \right]``

The function returns a real number if `O` is of `Hermitian` type or `Symmetric` type, and returns a complex number otherwise. You can make an operator `O` hermitian by using `Hermitian(O)`.

!!! note "List of observables and states"
    The observable `O` and state `ψ` can be given as a list of [`QuantumObject`](@ref), it returns a list of expectation values. If both of them are given as a list, it returns a `Matrix` of expectation values.

# Examples

```jldoctest
julia> ψ1 = 1 / √2 * (fock(10,2) + fock(10,4));

julia> ψ2 = coherent(10, 0.6 + 0.8im);

julia> a = destroy(10);

julia> expect(a' * a, ψ1) |> round
3.0 + 0.0im

julia> expect(Hermitian(a' * a), ψ1) |> round
3.0

julia> round.(expect([a' * a, a' + a, a], [ψ1, ψ2]), digits = 1)
3×2 Matrix{ComplexF64}:
 3.0+0.0im  1.0+0.0im
 0.0+0.0im  1.2-0.0im
 0.0+0.0im  0.6+0.8im
```
"""
expect(O::AbstractQuantumObject{Operator}, ψ::QuantumObject{Ket}) = dot(ψ.data, O.data, ψ.data)
expect(O::AbstractQuantumObject{Operator}, ψ::QuantumObject{Bra}) = expect(O, ψ')
expect(O::QuantumObject{Operator}, ρ::QuantumObject{Operator}) = tr(O * ρ)
expect(
    O::QuantumObject{Operator,DimsType,<:Union{<:Hermitian{TF},<:Symmetric{TR}}},
    ψ::QuantumObject{Ket},
) where {DimsType<:AbstractDimensions,TF<:Number,TR<:Real} = real(dot(ψ.data, O.data, ψ.data))
expect(
    O::QuantumObject{Operator,DimsType,<:Union{<:Hermitian{TF},<:Symmetric{TR}}},
    ψ::QuantumObject{Bra},
) where {DimsType<:AbstractDimensions,TF<:Number,TR<:Real} = real(expect(O, ψ'))
expect(
    O::QuantumObject{Operator,DimsType,<:Union{<:Hermitian{TF},<:Symmetric{TR}}},
    ρ::QuantumObject{Operator},
) where {DimsType<:AbstractDimensions,TF<:Number,TR<:Real} = real(tr(O * ρ))
expect(
    O::AbstractVector{<:AbstractQuantumObject{Operator,DimsType,<:Union{<:Hermitian{TF},<:Symmetric{TR}}}},
    ρ::QuantumObject,
) where {DimsType<:AbstractDimensions,TF<:Number,TR<:Real} = expect.(O, Ref(ρ))
function expect(O::AbstractVector{<:AbstractQuantumObject{Operator}}, ρ::QuantumObject)
    result = Vector{ComplexF64}(undef, length(O))
    result .= expect.(O, Ref(ρ))
    return result
end
expect(O::AbstractQuantumObject{Operator}, ρ::AbstractVector{<:QuantumObject}) = expect.(Ref(O), ρ)
function expect(
    O::AbstractVector{<:AbstractQuantumObject{Operator,DimsType,<:Union{<:Hermitian{TF},<:Symmetric{TR}}}},
    ρ::AbstractVector{<:QuantumObject},
) where {DimsType<:AbstractDimensions,TF<:Number,TR<:Real}
    N_ops = length(O)
    result = Matrix{Float64}(undef, N_ops, length(ρ))
    for i in 1:N_ops
        result[i, :] .= expect.(Ref(O[i]), ρ)
    end
    return result
end
function expect(O::AbstractVector{<:AbstractQuantumObject{Operator}}, ρ::AbstractVector{<:QuantumObject})
    N_ops = length(O)
    result = Matrix{ComplexF64}(undef, N_ops, length(ρ))
    for i in 1:N_ops
        result[i, :] .= expect.(Ref(O[i]), ρ)
    end
    return result
end

@doc raw"""
    variance(O::QuantumObject, ψ::Union{QuantumObject,Vector{QuantumObject}})

Variance of the [`Operator`](@ref) `O`: ``\langle\hat{O}^2\rangle - \langle\hat{O}\rangle^2``,

where ``\langle\hat{O}\rangle`` is the expectation value of `O` with the state `ψ` (see also [`expect`](@ref)), and the state `ψ` can be a [`Ket`](@ref), [`Bra`](@ref) or [`Operator`](@ref).

The function returns a real number if `O` is hermitian, and returns a complex number otherwise.

Note that `ψ` can also be given as a list of [`QuantumObject`](@ref), it returns a list of expectation values.
"""
variance(O::QuantumObject{Operator}, ψ::QuantumObject) = expect(O^2, ψ) - expect(O, ψ)^2
variance(O::QuantumObject{Operator}, ψ::Vector{<:QuantumObject}) = expect(O^2, ψ) .- expect(O, ψ) .^ 2

@doc raw"""
    to_dense(A::QuantumObject)

Converts a sparse QuantumObject to a dense QuantumObject.
"""
to_dense(A::QuantumObject) = QuantumObject(to_dense(A.data), A.type, A.dimensions)
to_dense(A::MT) where {MT<:AbstractSparseArray} = Array(A)
to_dense(A::MT) where {MT<:AbstractArray} = A

to_dense(::Type{T}, A::AbstractSparseArray) where {T<:Number} = Array{T}(A)
to_dense(::Type{T1}, A::AbstractArray{T2}) where {T1<:Number,T2<:Number} = Array{T1}(A)
to_dense(::Type{T}, A::AbstractArray{T}) where {T<:Number} = A

function to_dense(::Type{M}) where {M<:Union{Diagonal,SparseMatrixCSC}}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return Matrix{par[1]}
end

to_dense(::Type{M}) where {M<:AbstractMatrix} = M

@doc raw"""
    to_sparse(A::QuantumObject)

Converts a dense QuantumObject to a sparse QuantumObject.
"""
to_sparse(A::QuantumObject, tol::Real = 1e-10) = QuantumObject(to_sparse(A.data, tol), A.type, A.dimensions)
function to_sparse(A::MT, tol::Real = 1e-10) where {MT<:AbstractMatrix}
    idxs = findall(@. abs(A) > tol)
    row_indices = getindex.(idxs, 1)
    col_indices = getindex.(idxs, 2)
    vals = getindex(A, idxs)
    return sparse(row_indices, col_indices, vals, size(A)...)
end
function to_sparse(A::VT, tol::Real = 1e-10) where {VT<:AbstractVector}
    idxs = findall(@. abs(A) > tol)
    vals = getindex(A, idxs)
    return sparsevec(idxs, vals, length(A))
end

@doc raw"""
    kron(A::AbstractQuantumObject, B::AbstractQuantumObject, ...)
    tensor(A::AbstractQuantumObject, B::AbstractQuantumObject, ...)
    ⊗(A::AbstractQuantumObject, B::AbstractQuantumObject, ...)
    A ⊗ B

Returns the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) ``\hat{A} \otimes \hat{B} \otimes \cdots``.

!!! note
    `tensor` and `⊗` (where `⊗` can be typed by tab-completing `\otimes` in the REPL) are synonyms of `kron`.

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

julia> O = kron(a, a);

julia> size(a), size(O)
((20, 20), (400, 400))

julia> a.dims, O.dims
([20], [20, 20])
```
"""
function Base.kron(
    A::AbstractQuantumObject{OpType,<:Dimensions},
    B::AbstractQuantumObject{OpType,<:Dimensions},
) where {OpType<:Union{Ket,Bra,Operator}}
    QType = promote_op_type(A, B)
    _lazy_tensor_warning(A.data, B.data)
    return QType(kron(A.data, B.data), A.type, Dimensions((A.dimensions.to..., B.dimensions.to...)))
end

# if A and B are both Operator but either one of them has GeneralDimensions
for ADimType in (:Dimensions, :GeneralDimensions)
    for BDimType in (:Dimensions, :GeneralDimensions)
        if !(ADimType == BDimType == :Dimensions) # not for this case because it's already implemented
            @eval begin
                function Base.kron(
                    A::AbstractQuantumObject{Operator,<:$ADimType},
                    B::AbstractQuantumObject{Operator,<:$BDimType},
                )
                    QType = promote_op_type(A, B)
                    _lazy_tensor_warning(A.data, B.data)
                    return QType(
                        kron(A.data, B.data),
                        Operator(),
                        GeneralDimensions(
                            (get_dimensions_to(A)..., get_dimensions_to(B)...),
                            (get_dimensions_from(A)..., get_dimensions_from(B)...),
                        ),
                    )
                end
            end
        end
    end
end

# if A and B are different type (must return Operator with GeneralDimensions)
for AOpType in (:Ket, :Bra, :Operator)
    for BOpType in (:Ket, :Bra, :Operator)
        if (AOpType != BOpType)
            @eval begin
                function Base.kron(A::AbstractQuantumObject{$AOpType}, B::AbstractQuantumObject{$BOpType})
                    QType = promote_op_type(A, B)
                    _lazy_tensor_warning(A.data, B.data)
                    return QType(
                        kron(A.data, B.data),
                        Operator(),
                        GeneralDimensions(
                            (get_dimensions_to(A)..., get_dimensions_to(B)...),
                            (get_dimensions_from(A)..., get_dimensions_from(B)...),
                        ),
                    )
                end
            end
        end
    end
end

Base.kron(A::AbstractQuantumObject) = A
function Base.kron(A::Vector{<:AbstractQuantumObject})
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
    vector_to_operator(A::QuantumObject)

Convert a quantum object from vector ([`OperatorKet`](@ref)-type) to matrix ([`Operator`](@ref)-type)

!!! note
    `vector_to_operator` is a synonym of `vec2mat`.
"""
vec2mat(A::QuantumObject{OperatorKet}) = QuantumObject(vec2mat(A.data), Operator(), A.dimensions)

@doc raw"""
    mat2vec(A::QuantumObject)
    operator_to_vector(A::QuantumObject)

Convert a quantum object from matrix ([`Operator`](@ref)-type) to vector ([`OperatorKet`](@ref)-type)

!!! note
    `operator_to_vector` is a synonym of `mat2vec`.
"""
mat2vec(A::QuantumObject{Operator}) = QuantumObject(mat2vec(A.data), OperatorKet(), A.dimensions)

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
