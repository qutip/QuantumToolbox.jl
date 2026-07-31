#=
LinearAlgebra Eigen solvers and results for QuantumObject
=#

export EigsolveResult

@doc raw"""
    struct EigsolveResult

A struct containing the eigenvalues, the eigenvectors, and some information from the solver

# Fields (Attributes)
- `values::AbstractVector`: the eigenvalues
- `vectors::AbstractMatrix`: the transformation matrix (eigenvectors)
- `type::Union{Nothing,QuantumObjectType}`: the type of [`QuantumObject`](@ref), or `nothing` means solving eigen equation for general matrix
- `dimensions::Union{Nothing,Dimensions}`: the `dimensions` of [`QuantumObject`](@ref), or `nothing` means solving eigen equation for general matrix
- `iter::Int`: the number of iteration during the solving process
- `numops::Int` : number of times the linear map was applied in krylov methods
- `converged::Bool`: Whether the result is converged

!!! note "`dims` property"
    For a given `res::EigsolveResult`, `res.dims` or `getproperty(res, :dims)` returns its `dimensions` in the type of integer-vector.

# Examples
One can obtain the eigenvalues and the corresponding [`QuantumObject`](@ref)-type eigenvectors by:
```jldoctest
julia> result = eigen(sigmax())
EigsolveResult:   type=Operator()   dims=([2], [2])
values:
2-element Vector{Float64}:
 -1.0
  1.0
vectors:
2×2 Matrix{ComplexF64}:
 -0.707107+0.0im  0.707107+0.0im
  0.707107+0.0im  0.707107+0.0im

julia> λ, ψ, U = result;

julia> λ # eigenvalues
2-element Vector{Float64}:
 -1.0
  1.0

julia> ψ # eigenvectors
2-element Vector{QuantumObject{Ket, Dimensions{Space, Space}, Vector{ComplexF64}}}:

Quantum Object:   type=Ket()   dims=([2], [1])   size=(2,)
2-element Vector{ComplexF64}:
 -0.7071067811865475 + 0.0im
  0.7071067811865475 + 0.0im

Quantum Object:   type=Ket()   dims=([2], [1])   size=(2,)
2-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im

julia> U # the transformation matrix
2×2 Matrix{ComplexF64}:
 -0.707107+0.0im  0.707107+0.0im
  0.707107+0.0im  0.707107+0.0im
```
"""
struct EigsolveResult{
        T1 <: AbstractVector{<:Number},
        T2 <: AbstractMatrix{<:Number},
        ObjType <: Union{Nothing, Operator, SuperOperator},
        DimType <: Union{Nothing, Dimensions},
    }
    values::T1
    vectors::T2
    type::ObjType
    dimensions::DimType
    iter::Int
    numops::Int
    converged::Bool
end

function Base.getproperty(res::EigsolveResult, key::Symbol)
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(res, :dimensions))
    else
        return getfield(res, key)
    end
end

Base.iterate(res::EigsolveResult) = (res.values, Val(:vector_list))
Base.iterate(res::EigsolveResult{T1, T2, Nothing}, ::Val{:vector_list}) where {T1, T2} =
    ([res.vectors[:, k] for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1, T2, Operator}, ::Val{:vector_list}) where {T1, T2} =
    ([QuantumObject(res.vectors[:, k], Ket(), Dimensions(res.dimensions.to, Space(1))) for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1, T2, SuperOperator}, ::Val{:vector_list}) where {T1, T2} =
    ([QuantumObject(res.vectors[:, k], OperatorKet(), Dimensions(res.dimensions.to, Space(1))) for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult, ::Val{:vectors}) = (res.vectors, Val(:done))
Base.iterate(res::EigsolveResult, ::Val{:done}) = nothing

function Base.show(io::IO, res::EigsolveResult)
    println(io, "EigsolveResult:   type=", res.type, "   dims=", _get_dims_string(res.dimensions))
    println(io, "values:")
    show(io, MIME("text/plain"), res.values)
    print(io, "\n")
    println(io, "vectors:")
    return show(io, MIME("text/plain"), res.vectors)
end

@doc raw"""
    LinearAlgebra.eigen(A::QuantumObject; kwargs...)

Calculates the eigenvalues and eigenvectors of the [`QuantumObject`](@ref) `A` using
the Julia [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) package.

```jldoctest
julia> a = destroy(5);

julia> H = a + a';

julia> using LinearAlgebra;

julia> E, ψ = eigen(H); # eigenvalues and eigenvectors

julia> round.(E, digits = 5)
5-element Vector{Float64}:
 -2.85697
 -1.35563
  0.0
  1.35563
  2.85697

julia> expect(H, ψ[1]) ≈ E[1]
true
```
"""
function LinearAlgebra.eigen(A::QuantumObject{OpType}; kwargs...) where {OpType <: Union{Operator, SuperOperator}}
    # This creates a weak Union type on CPU. See https://github.com/JuliaLang/LinearAlgebra.jl/issues/1498
    F = eigen(to_dense(A.data); kwargs...)

    E = F.values
    U = F.vectors
    settings.auto_tidyup && tidyup!(U)

    return EigsolveResult(E, U, A.type, A.dimensions, 0, 0, true)
end

@doc raw"""
    LinearAlgebra.eigvals(A::QuantumObject; kwargs...)

Same as [`eigen(A::QuantumObject; kwargs...)`](@ref) but for only the eigenvalues.
"""
LinearAlgebra.eigvals(A::QuantumObject{OpType}; kwargs...) where {OpType <: Union{Operator, SuperOperator}} =
    eigvals(to_dense(A.data); kwargs...)
