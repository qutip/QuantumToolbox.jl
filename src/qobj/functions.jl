#=
Functions which manipulates QuantumObject
=#

export ket2dm
export expect, variance
export to_dense, to_sparse
export multisite_operator
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

julia> ψ2 = coherent(10, 0.6);

julia> a = destroy(10);

julia> expect(a' * a, ψ1) |> real |> round
3.0

julia> expect([a' * a, a' + a, a], [ψ1, ψ2]) |> real
3×2 Matrix{Float64}:
 3.0  0.36
 0.0  1.2
 0.0  0.6
```
"""
expect(O::AbstractQuantumObject{Operator}, ψ::QuantumObject{Ket}) = dot(ψ, O, ψ) # check_mul_dimensions in dot
expect(O::AbstractQuantumObject{Operator}, ψ::QuantumObject{Bra}) = expect(O, ψ')
expect(O::QuantumObject{Operator}, ρ::QuantumObject{Operator}) = tr(O * ρ) # check_mul_dimensions in :(*)
expect(
    O::QuantumObject{Operator, DimsType, <:Union{<:Hermitian{TF}, <:Symmetric{TR}}},
    ψ::QuantumObject{Ket},
) where {DimsType <: Dimensions, TF <: Number, TR <: Real} = real(dot(ψ, O, ψ)) # check_mul_dimensions in dot
expect(
    O::QuantumObject{Operator, DimsType, <:Union{<:Hermitian{TF}, <:Symmetric{TR}}},
    ψ::QuantumObject{Bra},
) where {DimsType <: Dimensions, TF <: Number, TR <: Real} = real(expect(O, ψ'))
expect(
    O::QuantumObject{Operator, DimsType, <:Union{<:Hermitian{TF}, <:Symmetric{TR}}},
    ρ::QuantumObject{Operator},
) where {DimsType <: Dimensions, TF <: Number, TR <: Real} = real(tr(O * ρ)) # check_mul_dimensions in :(*)
expect(
    O::AbstractVector{<:AbstractQuantumObject{Operator, DimsType, <:Union{<:Hermitian{TF}, <:Symmetric{TR}}}},
    ρ::QuantumObject,
) where {DimsType <: Dimensions, TF <: Number, TR <: Real} = expect.(O, Ref(ρ))
function expect(O::AbstractVector{<:AbstractQuantumObject{Operator}}, ρ::QuantumObject)
    result = Vector{ComplexF64}(undef, length(O))
    result .= expect.(O, Ref(ρ))
    return result
end
expect(O::AbstractQuantumObject{Operator}, ρ::AbstractVector{<:QuantumObject}) = expect.(Ref(O), ρ)
function expect(
        O::AbstractVector{<:AbstractQuantumObject{Operator, DimsType, <:Union{<:Hermitian{TF}, <:Symmetric{TR}}}},
        ρ::AbstractVector{<:QuantumObject},
    ) where {DimsType <: Dimensions, TF <: Number, TR <: Real}
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
to_dense(A::MT) where {MT <: AbstractSparseArray} = Array(A)
to_dense(A::MT) where {MT <: AbstractArray} = A
to_dense(A::Diagonal) = diagm(A.diag)

to_dense(::Type{T}, A::AbstractSparseArray) where {T <: Number} = Array{T}(A)
to_dense(::Type{T1}, A::AbstractArray{T2}) where {T1 <: Number, T2 <: Number} = Array{T1}(A)
to_dense(::Type{T}, A::AbstractArray{T}) where {T <: Number} = A
to_dense(::Type{T}, A::Diagonal{T}) where {T <: Number} = diagm(A.diag)

function to_dense(::Type{M}) where {M <: Union{Diagonal, SparseMatrixCSC}}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return Matrix{par[1]}
end

to_dense(::Type{M}) where {M <: AbstractMatrix} = M

@doc raw"""
    to_sparse(A::QuantumObject, tol::Real = settings.tidyup_tol)

Converts a dense QuantumObject to a sparse QuantumObject by removing elements with absolute values smaller than `tol`.
"""
to_sparse(A::QuantumObject, tol::Real = settings.tidyup_tol) = QuantumObject(to_sparse(A.data, tol), A.type, A.dimensions)
function to_sparse(A::MT, tol::Real = settings.tidyup_tol) where {MT <: AbstractMatrix}
    idxs = findall(@. abs(A) > tol)
    row_indices = getindex.(idxs, 1)
    col_indices = getindex.(idxs, 2)
    vals = getindex(A, idxs)
    return sparse(row_indices, col_indices, vals, size(A)...)
end
function to_sparse(A::VT, tol::Real = settings.tidyup_tol) where {VT <: AbstractVector}
    idxs = findall(@. abs(A) > tol)
    vals = getindex(A, idxs)
    return sparsevec(idxs, vals, length(A))
end

to_sparse_if_needed(::Val{needed}, A::Union{QuantumObject, AbstractArray}, tol::Real = settings.tidyup_tol) where {needed} =
    needed ? to_sparse(A, tol) : A

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

Quantum Object:   type=Operator()   dims=([20], [20])   size=(20, 20)   ishermitian=false
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
(([20], [20]), ([20, 20], [20, 20]))
```
"""
Base.kron(A::AbstractQuantumObject) = A

# kron for two quantum objects A and B
for AOpType in (:Ket, :Bra, :Operator)
    for BOpType in (:Ket, :Bra, :Operator)
        # handle the expressions for different type-cases here
        # so that we don't need to evaluate if-condition in the function body
        KronOpType = (AOpType == BOpType) ? AOpType : :Operator
        dimensions_to_expr = (AOpType == BOpType == :Bra) ? :(Space(1)) : :(kron(A.dimensions.to, B.dimensions.to))
        dimensions_from_expr = (AOpType == BOpType == :Ket) ? :(Space(1)) : :(kron(A.dimensions.from, B.dimensions.from))

        @eval begin
            function Base.kron(A::AbstractQuantumObject{$AOpType}, B::AbstractQuantumObject{$BOpType})
                QType = promote_op_type(A, B)
                _lazy_tensor_warning(A.data, B.data)
                return QType(
                    kron(A.data, B.data),
                    $(KronOpType)(),
                    Dimensions(
                        $dimensions_to_expr,
                        $dimensions_from_expr,
                    ),
                )
            end
        end
    end
end
function Base.kron(A::Vector{<:AbstractQuantumObject})
    @warn "`tensor(A)` or `kron(A)` with `A` is a `Vector` can hurt performance. Try to use `tensor(A...)` or `kron(A...)` instead."
    return kron(A...)
end

@doc raw"""
    multisite_operator(dims::Union{AbstractVecOrTuple,Integer,Val}, pairs::Pair{Integer,QuantumObject{Operator}}...)

A Julia function for generating a multi-site operator.

For example, a ``N``-site operator ``\hat{O}`` with operators ``\hat{A}``, ``\hat{B}``, and ``\hat{C}`` acting on sites ``i``, ``j``, and ``k``, respectively, can be expressed as

```math
\begin{aligned}
\hat{O} &= \hat{A}_i \hat{B}_j \hat{C}_k\\
&= \left(\bigotimes_{n=1}^{i-1}\hat{\mathbb{1}}_n \right) \otimes \hat{A}_i \otimes \left(\bigotimes_{n=i+1}^{j-1}\hat{\mathbb{1}}_n \right) \otimes \hat{B}_j \otimes \left(\bigotimes_{n=j+1}^{k-1}\hat{\mathbb{1}}_n \right) \otimes \hat{C}_k \otimes \left(\bigotimes_{n=k+1}^{N}\hat{\mathbb{1}}_n \right).
\end{aligned}
```

# Arguments
- `dims::Union{AbstractVecOrTuple,Integer,Val}`: A list of integers representing the Hilbert space dimensions of each site. If `dims` is specified as an `Integer` or `Val`, it is assumed that all sites have the same Hilbert space dimension.
- `pairs::Pair{Integer,QuantumObject{Operator}}...`: A list of pairs where the first element of the pair is the site index and the second element is the [`Operator`](@ref) acting on that site.

# Returns
`QuantumObject`: A `QuantumObject` representing the multi-site operator.

# Examples

Consider two qubits ``\mathrm{A}`` and ``\mathrm{B}`` together with two cavity modes ``\mathrm{C1}`` and ``\mathrm{C2}`` (set both Hilbert space cutoffs to `5`). The following operators

- ``\hat{H}_\mathrm{AB} = \hat{\sigma}^x_\mathrm{A} \hat{\sigma}^x_\mathrm{B}``
- ``\hat{H}_\mathrm{AC1} = \hat{\sigma}^-_\mathrm{A} \hat{a}^\dagger_\mathrm{C1} + \hat{\sigma}^+_\mathrm{A} \hat{a}_\mathrm{C1}``
- ``\hat{H}_\mathrm{BC2} = \hat{\sigma}^-_\mathrm{B} \hat{a}^\dagger_\mathrm{C2} + \hat{\sigma}^+_\mathrm{B} \hat{a}_\mathrm{C2}``

can be generated by:

```jldoctest
julia> dims = (5, 2, 2, 5); # the Hilbert space of (C1, A, B, C2)

julia> σx = sigmax();

julia> σm = sigmam();

julia> a  = destroy(5);

julia> H_AB  = multisite_operator(dims, 2=>σx, 3=>σx);

julia> H_AC1 = multisite_operator(dims, 2=>σm, 1=>a') + multisite_operator(dims, 2=>σm', 1=>a);

julia> H_BC2 = multisite_operator(dims, 3=>σm, 4=>a') + multisite_operator(dims, 3=>σm', 4=>a);

julia> H_AB.dims == H_AC1.dims == H_BC2.dims == ([5, 2, 2, 5], [5, 2, 2, 5])
true
```

If all sites have the same Hilbert space dimension, you can just specify the number of sites:

```jldoctest
julia> N = Val(4); # four sites

julia> a = multisite_operator(N, 2=>destroy(3));

julia> a.dims
([3, 3, 3, 3], [3, 3, 3, 3])
```

Another example, for an `8`-site spin-``\frac{1}{2}`` chain, the operator ``\hat{\sigma}^x_5 \hat{\sigma}^z_7`` can be generated by

```jldoctest
julia> op = multisite_operator(Val(8), 5=>sigmax(), 7=>sigmaz());

julia> op.dims
([2, 2, 2, 2, 2, 2, 2, 2], [2, 2, 2, 2, 2, 2, 2, 2])
```
"""
function multisite_operator(dims::AbstractVecOrTuple{T}, pairs::Pair{<:Integer, <:QuantumObject{Operator}}...) where {T <: Integer}
    isempty(pairs) && throw(ArgumentError("At least one Pair of `site-index => operator` must be provided."))

    N = length(dims) # total number of sites
    sites_unsorted = collect(getfield.(pairs, :first))
    all(i -> 1 <= i <= N, sites_unsorted) || throw(ArgumentError("There are totally $N-sites, so site indices must satisfy 1 ≤ i ≤ $N."))

    idxs = sortperm(sites_unsorted)
    _sites = sites_unsorted[idxs]
    _ops = collect(getfield.(pairs, :second))[idxs]
    _dims = collect(dims) # Use this instead of a Tuple, to avoid type instability when indexing on a slice

    sites, ops, ElType = _get_unique_sites_ops_type(_sites, _ops)

    (all(isendomorphic, ops) && (_dims[sites] == [get_size(op.dimensions)[1] for op in ops])) ||
        throw(ArgumentError("The dimensions of the operators do not match the site dimensions."))

    data = kron(Eye{ElType}(prod(_dims[1:(sites[1] - 1)])), ops[1].data)
    for i in 2:length(sites)
        data = kron(data, Eye{ElType}(prod(_dims[(sites[i - 1] + 1):(sites[i] - 1)])), ops[i].data)
    end
    data = kron(data, Eye{ElType}(prod(_dims[(sites[end] + 1):end])))

    return QuantumObject(data; type = Operator(), dims = dims)
end
function multisite_operator(N::Union{Integer, Val}, pairs::Pair{<:Integer, <:QuantumObject{Operator}}...)
    isempty(pairs) && throw(ArgumentError("At least one Pair of `site-index => operator` must be provided."))

    d = get_size(first(pairs).second.dimensions)[1]
    dims = ntuple(j -> d, makeVal(N))
    return multisite_operator(dims, pairs...)
end

function _get_unique_sites_ops_type(sites, ops)
    unique_sites = unique(sites)
    unique_ops = map(i -> prod(ops[findall(==(i), sites)]), unique_sites)
    T = mapreduce(eltype, promote_type, unique_ops)

    return unique_sites, unique_ops, T
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
function vec2mat(A::QuantumObject{OperatorKet, <:Dimensions{<:LiouvilleSpace, Space}})
    op_dims = A.dimensions.to.op_dims
    m, n = get_size(op_dims)
    return QuantumObject(reshape(A.data, m, n), Operator(), op_dims)
end

@doc raw"""
    mat2vec(A::QuantumObject)
    operator_to_vector(A::QuantumObject)

Convert a quantum object from matrix ([`Operator`](@ref)-type) to vector ([`OperatorKet`](@ref)-type)

!!! note
    `operator_to_vector` is a synonym of `mat2vec`.
"""
mat2vec(A::QuantumObject{Operator}) = QuantumObject(mat2vec(A.data), OperatorKet(), Dimensions(LiouvilleSpace(A.dimensions), Space(1)))

@doc raw"""
    mat2vec(A::AbstractMatrix)

Converts a matrix to a vector.
"""
mat2vec(A::MT) where {MT <: AbstractMatrix} = vec(A) # reshape(A, :)
function mat2vec(A::MT) where {MT <: AbstractSparseMatrix}
    i, j, v = findnz(A)
    return sparsevec(i .+ (j .- 1) .* size(A, 1), v, prod(size(A)))
end
for op in (:Transpose, :Adjoint)
    @eval mat2vec(A::$op{T, <:AbstractSparseMatrix}) where {T <: Number} = mat2vec(sparse(A))
    @eval mat2vec(A::$op{T, MT}) where {T <: Number, MT <: AbstractMatrix} = mat2vec(MT(A))
end

function mat2vec(::Type{M}) where {M <: DenseMatrix}
    T = hasproperty(M, :body) ? M.body : M
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $M is not supported.")
    if npar == 2
        S = T.name.wrapper{par[1], 1}
    else
        S = T.name.wrapper{par[1], 1, par[3]}
    end
    return S
end

function mat2vec(::Type{M}) where {M <: SparseMatrixCSC}
    T = M
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return SparseVector{par[1], par[2]}
end

function mat2vec(::Type{M}) where {M <: Union{Adjoint{<:Number, <:SparseMatrixCSC}, Transpose{<:Number, <:SparseMatrixCSC}}}
    T = M.parameters[2]
    par = T.parameters
    npar = length(par)
    (2 == npar) || error("Type $M is not supported.")
    return SparseVector{par[1], par[2]}
end
