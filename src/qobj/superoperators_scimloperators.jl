#=
Custom AbstractSciMLOperator types for non-vectorized (matrix-form) superoperator actions.
These operate on density matrices in matrix form (N×N) rather than vectorized (N²) form.
=#

export SpostSuperOperator, SprePostSuperOperator

@doc raw"""
    SpostSuperOperator{T, OpType <: Union{AbstractMatrix, AbstractSciMLOperator}} <: AbstractSciMLOperator{T}

Represents the superoperator ``\mathcal{O}(\hat{B})[\hat{\rho}] = \hat{\rho}\hat{B}`` acting on the **non-vectorized**
density operator matrix ``\hat{\rho}``.

The action is defined as:
```math
\mathcal{O}(\hat{B})[\hat{\rho}] = \hat{\rho} \hat{B}
```

# Fields
- `R::OpType`: The operator ``\hat{B}`` (right-multiplier) as an `AbstractSciMLOperator`

See also [`SprePostSuperOperator`](@ref), [`spost`](@ref).
"""
struct SpostSuperOperator{T, OpType <: Union{AbstractMatrix, AbstractSciMLOperator}} <: AbstractSciMLOperator{T}
    R::OpType
    function SpostSuperOperator(R::OpType) where {OpType <: Union{AbstractMatrix, AbstractSciMLOperator}}
        T = eltype(R)
        return new{T, OpType}(R)
    end
end

function Base.show(io::IO, L::SpostSuperOperator)
    a, b = size(L)
    return print(io, "SpostSuperOperator($a × $b)")
end

Base.size(op::SpostSuperOperator) = size(op.R)

function LinearAlgebra.mul!(v::AbstractMatrix, op::SpostSuperOperator, u::AbstractMatrix)
    mul!(v, u, op.R)
    return v
end

function LinearAlgebra.mul!(v::AbstractMatrix, op::SpostSuperOperator, u::AbstractMatrix, α::Number, β::Number)
    mul!(v, u, op.R, α, β)
    return v
end

function Base.:*(A::SpostSuperOperator, B::AbstractMatrix)
    C = similar(B, size(A, 1), size(B, 2))
    mul!(C, A, B)
    return C
end

SciMLOperators.isconstant(op::SpostSuperOperator) = isconstant(op.R)
SciMLOperators.iscached(op::SpostSuperOperator) = true  # no cache needed

function SciMLOperators.cache_operator(op::SpostSuperOperator, u::AbstractMatrix)
    R_cached = cache_operator(op.R, u)
    return SpostSuperOperator(R_cached)
end

function SciMLOperators.update_coefficients!(op::SpostSuperOperator, u::AbstractMatrix, p, t)
    update_coefficients!(op.R, u, p, t)
    return op
end

@doc raw"""
    SprePostSuperOperator{T, LOp <: Union{AbstractMatrix, AbstractSciMLOperator}, ROp <: Union{AbstractMatrix, AbstractSciMLOperator}, CT} <: AbstractSciMLOperator{T}

Represents the superoperator ``\mathcal{O}(\hat{A}, \hat{B})[\hat{\rho}] = \hat{A}\hat{\rho}\hat{B}`` acting on the
**non-vectorized** density operator matrix ``\hat{\rho}``.

The action is defined as:
```math
\mathcal{O}(\hat{A}, \hat{B})[\hat{\rho}] = \hat{A} \hat{\rho} \hat{B}
```

# Fields
- `L::LOp`: The operator ``\hat{A}`` (left-multiplier) as an `AbstractSciMLOperator`
- `R::ROp`: The operator ``\hat{B}`` (right-multiplier) as an `AbstractSciMLOperator`
- `cache::CT`: A cache matrix for intermediate computation (allocated via [`cache_operator`](@ref))

See also [`SpostSuperOperator`](@ref), [`sprepost`](@ref).
"""
struct SprePostSuperOperator{T, LOp <: Union{AbstractMatrix, AbstractSciMLOperator}, ROp <: Union{AbstractMatrix, AbstractSciMLOperator}, CT} <: AbstractSciMLOperator{T}
    L::LOp
    R::ROp
    cache::CT
    function SprePostSuperOperator(L::LOp, R::ROp, cache::CT = nothing) where {LOp <: Union{AbstractMatrix, AbstractSciMLOperator}, ROp <: Union{AbstractMatrix, AbstractSciMLOperator}, CT}
        T = promote_type(eltype(L), eltype(R))
        return new{T, LOp, ROp, CT}(L, R, cache)
    end
end

function Base.show(io::IO, L::SprePostSuperOperator)
    a, b = size(L)
    return print(io, "SprePostSuperOperator($a × $b)")
end

Base.size(op::SprePostSuperOperator) = size(op.L)

function LinearAlgebra.mul!(v::AbstractMatrix, op::SprePostSuperOperator, u::AbstractMatrix)
    iscached(op) || throw(ArgumentError("The cache for the SprePostSuperOperator must be initialized before multiplication. Use `cache_operator` to initialize the cache."))
    mul!(op.cache, op.L, u)     # cache = L * u
    mul!(v, op.cache, op.R)     # v = cache * R = L * u * R
    return v
end

function LinearAlgebra.mul!(v::AbstractMatrix, op::SprePostSuperOperator, u::AbstractMatrix, α::Number, β::Number)
    iscached(op) || throw(ArgumentError("The cache for the SprePostSuperOperator must be initialized before multiplication. Use `cache_operator` to initialize the cache."))
    mul!(op.cache, op.L, u)            # cache = L * u
    mul!(v, op.cache, op.R, α, β)     # v = α * (L * u * R) + β * v
    return v
end

function Base.:*(A::SprePostSuperOperator, B::AbstractMatrix)
    iscached(A) || throw(ArgumentError("The cache for the SprePostSuperOperator must be initialized before multiplication. Use `cache_operator` to initialize the cache."))
    C = similar(B, size(A, 1), size(B, 2))
    mul!(C, A, B)
    return C
end

SciMLOperators.isconstant(op::SprePostSuperOperator) = isconstant(op.L) && isconstant(op.R)
SciMLOperators.iscached(op::SprePostSuperOperator) = op.cache !== nothing

function SciMLOperators.update_coefficients!(op::SprePostSuperOperator, u::AbstractMatrix, p, t)
    update_coefficients!(op.L, u, p, t)
    update_coefficients!(op.R, u, p, t)
    return op
end

function SciMLOperators.cache_operator(op::SprePostSuperOperator, u::AbstractMatrix)
    L_cached = cache_operator(op.L, u)
    R_cached = cache_operator(op.R, u)
    cache = similar(u)
    return SprePostSuperOperator(L_cached, R_cached, cache)
end

SciMLOperators.getcache(op::SprePostSuperOperator) = op.cache

function SciMLOperators._get_cache_shapes(::SprePostSuperOperator, u::AbstractMatrix)
    return size(u)
end
