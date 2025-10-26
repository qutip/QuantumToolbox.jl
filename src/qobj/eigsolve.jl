#=
Eigen solvers and results for QuantumObject
=#

export EigsolveResult
export eigenenergies, eigenstates, eigsolve
export eigsolve_al

@doc raw"""
    struct EigsolveResult

A struct containing the eigenvalues, the eigenvectors, and some information from the solver

# Fields (Attributes)
- `values::AbstractVector`: the eigenvalues
- `vectors::AbstractMatrix`: the transformation matrix (eigenvectors)
- `type::Union{Nothing,QuantumObjectType}`: the type of [`QuantumObject`](@ref), or `nothing` means solving eigen equation for general matrix
- `dimensions::Union{Nothing,AbstractDimensions}`: the `dimensions` of [`QuantumObject`](@ref), or `nothing` means solving eigen equation for general matrix
- `iter::Int`: the number of iteration during the solving process
- `numops::Int` : number of times the linear map was applied in krylov methods
- `converged::Bool`: Whether the result is converged

!!! note "`dims` property"
    For a given `res::EigsolveResult`, `res.dims` or `getproperty(res, :dims)` returns its `dimensions` in the type of integer-vector.

# Examples
One can obtain the eigenvalues and the corresponding [`QuantumObject`](@ref)-type eigenvectors by:
```jldoctest
julia> result = eigenstates(sigmax())
EigsolveResult:   type=Operator()   dims=[2]
values:
2-element Vector{ComplexF64}:
 -1.0 + 0.0im
  1.0 + 0.0im
vectors:
2×2 Matrix{ComplexF64}:
 -0.707107+0.0im  0.707107+0.0im
  0.707107+0.0im  0.707107+0.0im

julia> λ, ψ, U = result;

julia> λ
2-element Vector{ComplexF64}:
 -1.0 + 0.0im
  1.0 + 0.0im

julia> ψ
2-element Vector{QuantumObject{Ket, Dimensions{1, Tuple{Space}}, Vector{ComplexF64}}}:

Quantum Object:   type=Ket()   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 -0.7071067811865475 + 0.0im
  0.7071067811865475 + 0.0im

Quantum Object:   type=Ket()   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im

julia> U
2×2 Matrix{ComplexF64}:
 -0.707107+0.0im  0.707107+0.0im
  0.707107+0.0im  0.707107+0.0im
```
"""
struct EigsolveResult{
    T1<:Vector{<:Number},
    T2<:AbstractMatrix{<:Number},
    ObjType<:Union{Nothing,Operator,SuperOperator},
    DimType<:Union{Nothing,AbstractDimensions},
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
Base.iterate(res::EigsolveResult{T1,T2,Nothing}, ::Val{:vector_list}) where {T1,T2} =
    ([res.vectors[:, k] for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1,T2,Operator}, ::Val{:vector_list}) where {T1,T2} =
    ([QuantumObject(res.vectors[:, k], Ket(), res.dimensions) for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1,T2,SuperOperator}, ::Val{:vector_list}) where {T1,T2} =
    ([QuantumObject(res.vectors[:, k], OperatorKet(), res.dimensions) for k in 1:length(res.values)], Val(:vectors))
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

struct ArnoldiLindbladIntegratorMap{T,TS,TI} <: AbstractLinearMap{T,TS}
    elty::Type{T}
    size::TS
    integrator::TI
end

function LinearAlgebra.mul!(y::AbstractVector, A::ArnoldiLindbladIntegratorMap, x::AbstractVector)
    reinit!(A.integrator, x)
    solve!(A.integrator)
    return copyto!(y, A.integrator.u)
end

struct EigsolveInverseMap{T,TS,TI} <: AbstractLinearMap{T,TS}
    elty::Type{T}
    size::TS
    linsolve::TI
end

function LinearAlgebra.mul!(y::AbstractVector, A::EigsolveInverseMap, x::AbstractVector)
    A.linsolve.b .= x
    return copyto!(y, solve!(A.linsolve).u)
end

function _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, m, β, sorted_vals, sortby, rev)
    copyto!(Uₘ, Hₘ)
    F = schur!(Uₘ)

    values = F.values
    sortperm!(sorted_vals, values, by = sortby, rev = rev)

    select = fill(false, length(values))
    @inbounds for j in 1:k
        select[sorted_vals[j]] = true
    end

    ordschur!(F, select)

    copyto!(Hₘ, F.T)
    copyto!(Uₘ, F.Z)
    mul!(f, Uₘᵥ, β)

    return nothing
end

# Pure Julia implementation of computing right eigenvectors from Schur form
# Instead of using LAPACK.trevc!('R', 'A', select, Tₘ)
function _schur_right_eigenvectors(Tₘ, k)
    n = size(Tₘ, 1)
    vecs = zeros(eltype(Tₘ), n, k)
    k == 0 && return vecs

    value_tol(λ) = eps(typeof(abs(λ))) * max(one(typeof(abs(λ))), abs(λ))

    @inbounds for col in 1:k
        vec = view(vecs, :, col)
        fill!(vec, zero(eltype(Tₘ)))
        vec[col] = one(eltype(Tₘ))
        λ = Tₘ[col, col]

        for row in (col-1):-1:1
            acc = zero(eltype(Tₘ))
            for inner in (row + 1):col
                acc += Tₘ[row, inner] * vec[inner]
            end
            denom = Tₘ[row, row] - λ
            if abs(denom) <= value_tol(λ)
                vec[row] = zero(eltype(Tₘ))
            else
                vec[row] = -acc / denom
            end
        end

        LinearAlgebra.normalize!(vec)
    end

    return vecs
end

function _eigsolve(
    A,
    b::AbstractVector{T},
    type::ObjType,
    dimensions::Union{Nothing,AbstractDimensions},
    k::Int = 1,
    m::Int = max(20, 2 * k + 1);
    tol::Real = 1e-8,
    maxiter::Int = 200,
    sortby::Function = abs2,
    rev = true,
) where {T<:Number,ObjType<:Union{Nothing,Operator,SuperOperator}}
    n = size(A, 2)
    V = similar(b, n, m + 1)
    H = zeros(T, m + 1, m)

    arnoldi_init!(A, b, V, H)

    for i in 2:m
        β = arnoldi_step!(A, V, H, i)
        if β < tol && i > k
            m = i # happy breakdown
            break
        end
    end

    f = ones(eltype(A), m)

    Vₘ = view(V, :, 1:m)
    Hₘ = view(H, 1:m, 1:m)
    qₘ = view(V, :, m + 1)
    βeₘ = view(H, m + 1, 1:m)
    β = real(H[m+1, m])
    Uₘ = one(Hₘ)

    Uₘᵥ = view(Uₘ, m, 1:m)

    cache0 = similar(b, m, m)
    cache1 = similar(b, size(V, 1), m)
    cache2 = similar(H, m)
    sorted_vals = Array{Int}(undef, m)

    V₁ₖ = view(V, :, 1:k)
    Vₖ₊₁ = view(V, :, k + 1)
    Hₖ₊₁₁ₖ = view(H, k + 1, 1:k)
    cache0₁ₖ = view(cache0, :, 1:k)
    cache1₁ₖ = view(cache1, :, 1:k)
    cache2₁ₖ = view(cache2, 1:k)

    M = typeof(cache0)

    _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, m, β, sorted_vals, sortby, rev)

    numops = m
    iter = 0
    while iter < maxiter && count(x -> abs(x) < tol, f) < k && β > tol
        # println( A * Vₘ * Uₘ ≈ Vₘ * Uₘ * M(Tₘ) + qₘ * M(transpose(βeₘ)) * Uₘ )     # SHOULD BE TRUE

        copyto!(cache0, Uₘ)
        mul!(cache1, Vₘ, cache0)
        copyto!(V₁ₖ, cache1₁ₖ)
        copyto!(Vₖ₊₁, qₘ)
        mul!(cache2, transpose(Uₘ), βeₘ) # transpose(βeₘ) * Uₘ
        copyto!(Hₖ₊₁₁ₖ, cache2₁ₖ)

        # println( A * view(V, :, 1:k) ≈ view(V, :, 1:k) * M(view(H, 1:k, 1:k)) + qₘ * M(transpose(view(transpose(βeₘ) * Uₘ, 1:k))) )     # SHOULD BE TRUE

        for j in (k+1):m
            β = arnoldi_step!(A, V, H, j)
            if β < tol
                numops += j - k - 1
                break
            end
        end

        # println( A * Vₘ ≈ Vₘ * M(Hₘ) + qₘ * M(transpose(βeₘ)) )     # SHOULD BE TRUE

        _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, m, β, sorted_vals, sortby, rev)

        numops += m - k - 1
        iter += 1
    end

    Tₘ = Hₘ
    vals = diag(view(Tₘ, 1:k, 1:k))
    VR = _schur_right_eigenvectors(Tₘ, k)
    mul!(cache0₁ₖ, Uₘ, VR)
    mul!(cache1₁ₖ, Vₘ, cache0₁ₖ)
    vecs = copy(cache1₁ₖ)
    settings.auto_tidyup && tidyup!(vecs)

    return EigsolveResult(vals, vecs, type, dimensions, iter, numops, (iter < maxiter))
end

@doc raw"""
    eigsolve(A::QuantumObject; 
        v0::Union{Nothing,AbstractVector}=nothing, 
        sigma::Union{Nothing, Real}=nothing,
        eigvals::Int = 1,
        krylovdim::Int = max(20, 2*k+1),
        tol::Real = 1e-8,
        maxiter::Int = 200,
        solver::Union{Nothing, SciMLLinearSolveAlgorithm} = nothing,
        sortby::Function = abs2,
        rev::Bool = true,
        kwargs...)

Solve for the eigenvalues and eigenvectors of a matrix `A` using the Arnoldi method.

# Arguments
- `A::QuantumObject`: the [`QuantumObject`](@ref) to solve eigenvalues and eigenvectors.
- `v0::Union{Nothing,AbstractVector}`: the initial vector for the Arnoldi method. Default is a random vector.
- `sigma::Union{Nothing, Real}`: the shift for the eigenvalue problem. Default is `nothing`.
- `eigvals::Int`: the number of eigenvalues to compute. Default is `1`.
- `krylovdim::Int`: the dimension of the Krylov subspace. Default is `max(20, 2*k+1)`.
- `tol::Real`: the tolerance for the Arnoldi method. Default is `1e-8`.
- `maxiter::Int`: the maximum number of iterations for the Arnoldi method. Default is `200`.
- `solver::Union{Nothing, SciMLLinearSolveAlgorithm}`: the linear solver algorithm. Default is `nothing`.
- `sortby::Function`: the function to sort eigenvalues. Default is `abs2`.
- `rev::Bool`: whether to sort in descending order. Default is `true`.
- `kwargs`: Additional keyword arguments passed to the solver.

# Notes
- For more details about `solver` and extra `kwargs`, please refer to [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/)

# Returns
- `EigsolveResult`: A struct containing the eigenvalues, the eigenvectors, and some information about the eigsolver
"""
function eigsolve(
    A::QuantumObject;
    v0::Union{Nothing,AbstractVector} = nothing,
    sigma::Union{Nothing,Real} = nothing,
    eigvals::Int = 1,
    krylovdim::Int = max(20, 2 * eigvals + 1),
    tol::Real = 1e-8,
    maxiter::Int = 200,
    solver::Union{Nothing,SciMLLinearSolveAlgorithm} = nothing,
    sortby::Function = abs2,
    rev::Bool = true,
    kwargs...,
)
    return eigsolve(
        A.data;
        v0 = v0,
        type = A.type,
        dimensions = A.dimensions,
        sigma = sigma,
        eigvals = eigvals,
        krylovdim = krylovdim,
        tol = tol,
        maxiter = maxiter,
        solver = solver,
        sortby = sortby,
        rev = rev,
        kwargs...,
    )
end

function eigsolve(
    A;
    v0::Union{Nothing,AbstractVector} = nothing,
    type::Union{Nothing,Operator,SuperOperator} = nothing,
    dimensions = nothing,
    sigma::Union{Nothing,Real} = nothing,
    eigvals::Int = 1,
    krylovdim::Int = max(20, 2 * eigvals + 1),
    tol::Real = 1e-8,
    maxiter::Int = 200,
    solver::Union{Nothing,SciMLLinearSolveAlgorithm} = nothing,
    sortby::Function = abs2,
    rev::Bool = true,
    kwargs...,
)
    T = eltype(A)
    isH = ishermitian(A)
    v0 === nothing && (v0 = normalize!(rand(T, size(A, 1))))

    if sigma === nothing
        res = _eigsolve(
            A,
            v0,
            type,
            dimensions,
            eigvals,
            krylovdim,
            tol = tol,
            maxiter = maxiter,
            sortby = sortby,
            rev = rev,
        )
        vals = res.values
    else
        Aₛ = A - sigma * I
        solver === nothing && (solver = isH ? KrylovJL_MINRES() : KrylovJL_GMRES())

        kwargs2 = (; kwargs...)
        condition = !haskey(kwargs2, :Pl) && typeof(A) <: SparseMatrixCSC
        condition && (kwargs2 = merge(kwargs2, (Pl = ilu(Aₛ, τ = 0.01),)))

        !haskey(kwargs2, :abstol) && (kwargs2 = merge(kwargs2, (abstol = tol * 1e-6,)))
        !haskey(kwargs2, :reltol) && (kwargs2 = merge(kwargs2, (reltol = tol * 1e-6,)))

        prob = LinearProblem(Aₛ, v0)
        linsolve = init(prob, solver; kwargs2...)

        Amap = EigsolveInverseMap(T, size(A), linsolve)

        res = _eigsolve(
            Amap,
            v0,
            type,
            dimensions,
            eigvals,
            krylovdim,
            tol = tol,
            maxiter = maxiter,
            sortby = sortby,
            rev = rev,
        )
        vals = @. (1 + sigma * res.values) / res.values
    end

    vecs = res.vectors
    settings.auto_tidyup && tidyup!(vecs)

    return EigsolveResult(vals, vecs, res.type, res.dimensions, res.iter, res.numops, res.converged)
end

@doc raw"""
    eigsolve_al(
        H::Union{AbstractQuantumObject{HOpType},Tuple},
        T::Real,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        params::NamedTuple = NamedTuple(),
        ρ0::AbstractMatrix = rand_dm(prod(H.dimensions)).data,
        eigvals::Int = 1,
        krylovdim::Int = min(10, size(H, 1)),
        maxiter::Int = 200,
        eigstol::Real = 1e-6,
        sortby::Function = abs2,
        rev::Bool = true,
        kwargs...,
    )

Solve the eigenvalue problem for a Liouvillian superoperator `L` using the Arnoldi-Lindblad method.

# Arguments
- `H`: The Hamiltonian (or directly the Liouvillian) of the system. It can be a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a tuple of the form supported by [`mesolve`](@ref).
- `T`: The time at which to evaluate the time evolution.
- `c_ops`: A vector of collapse operators. Default is `nothing` meaning the system is closed.
- `alg`: The differential equation solver algorithm. Default is `Tsit5()`.
- `params`: A `NamedTuple` containing the parameters of the system.
- `ρ0`: The initial density matrix. If not specified, a random density matrix is used.
- `eigvals`: The number of eigenvalues to compute.
- `krylovdim`: The dimension of the Krylov subspace.
- `maxiter`: The maximum number of iterations for the eigsolver.
- `eigstol`: The tolerance for the eigsolver.
- `sortby::Function`: the function to sort eigenvalues. Default is `abs2`.
- `rev::Bool`: whether to sort in descending order. Default is `true`.
- `kwargs`: Additional keyword arguments passed to the differential equation solver.

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- `EigsolveResult`: A struct containing the eigenvalues, the eigenvectors, and some information about the eigsolver

# References
- [1] Minganti, F., & Huybrechts, D. (2022). Arnoldi-Lindblad time evolution: Faster-than-the-clock algorithm for the spectrum of time-independent and Floquet open quantum systems. Quantum, 6, 649.
"""
function eigsolve_al(
    H::Union{AbstractQuantumObject{HOpType},Tuple},
    T::Real,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    params::NamedTuple = NamedTuple(),
    ρ0::AbstractMatrix = rand_dm(prod(H.dimensions)).data,
    eigvals::Int = 1,
    krylovdim::Int = min(10, size(H, 1)),
    maxiter::Int = 200,
    eigstol::Real = 1e-6,
    sortby::Function = abs2,
    rev::Bool = true,
    kwargs...,
) where {HOpType<:Union{Operator,SuperOperator}}
    L_evo = _mesolve_make_L_QobjEvo(H, c_ops)
    prob = mesolveProblem(
        L_evo,
        QuantumObject(ρ0, type = Operator(), dims = H.dimensions),
        [zero(T), T];
        params = params,
        progress_bar = Val(false),
        kwargs...,
    ).prob
    integrator = init(prob, alg)

    Lmap = ArnoldiLindbladIntegratorMap(eltype(H), size(L_evo), integrator)

    res = _eigsolve(
        Lmap,
        mat2vec(ρ0),
        L_evo.type,
        L_evo.dimensions,
        eigvals,
        krylovdim,
        maxiter = maxiter,
        tol = eigstol,
        sortby = sortby,
        rev = rev,
    )

    vals = similar(res.values)
    vecs = similar(res.vectors)

    for i in eachindex(res.values)
        vec = view(res.vectors, :, i)
        vals[i] = dot(vec, L_evo.data, vec)
        @. vecs[:, i] = vec * exp(-1im * angle(vec[1]))
    end

    settings.auto_tidyup && tidyup!(vecs)

    return EigsolveResult(vals, vecs, res.type, res.dimensions, res.iter, res.numops, res.converged)
end

@doc raw"""
    LinearAlgebra.eigen(A::QuantumObject; kwargs...)

Calculates the eigenvalues and eigenvectors of the [`QuantumObject`](@ref) `A` using
the Julia [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) package.

```jldoctest
julia> a = destroy(5);

julia> H = a + a';

julia> using LinearAlgebra;

julia> E, ψ, U = eigen(H)
EigsolveResult:   type=Operator()   dims=[5]
values:
5-element Vector{ComplexF64}:
       -2.8569700138728 + 0.0im
    -1.3556261799742608 + 0.0im
 1.3322676295501878e-15 + 0.0im
     1.3556261799742677 + 0.0im
     2.8569700138728056 + 0.0im
vectors:
5×5 Matrix{ComplexF64}:
  0.106101+0.0im  -0.471249-0.0im  …   0.471249+0.0im  0.106101+0.0im
 -0.303127-0.0im   0.638838+0.0im      0.638838+0.0im  0.303127+0.0im
  0.537348+0.0im  -0.279149-0.0im      0.279149+0.0im  0.537348+0.0im
 -0.638838-0.0im  -0.303127-0.0im     -0.303127-0.0im  0.638838+0.0im
  0.447214+0.0im   0.447214+0.0im     -0.447214-0.0im  0.447214+0.0im

julia> expect(H, ψ[1]) ≈ E[1]
true
```
"""
function LinearAlgebra.eigen(A::QuantumObject{OpType}; kwargs...) where {OpType<:Union{Operator,SuperOperator}}
    MT = typeof(A.data)
    F = eigen(to_dense(A.data); kwargs...)
    # This fixes a type inference issue. But doesn't work for GPU arrays
    E::mat2vec(to_dense(MT)) = F.values
    U::to_dense(MT) = F.vectors
    settings.auto_tidyup && tidyup!(U)

    return EigsolveResult(E, U, A.type, A.dimensions, 0, 0, true)
end

@doc raw"""
    LinearAlgebra.eigvals(A::QuantumObject; kwargs...)

Same as [`eigen(A::QuantumObject; kwargs...)`](@ref) but for only the eigenvalues.
"""
LinearAlgebra.eigvals(A::QuantumObject{OpType}; kwargs...) where {OpType<:Union{Operator,SuperOperator}} =
    eigvals(to_dense(A.data); kwargs...)

@doc raw"""
    eigenenergies(A::QuantumObject; sparse::Bool=false, kwargs...)

Calculate the eigenenergies

# Arguments
- `A::QuantumObject`: the [`QuantumObject`](@ref) to solve eigenvalues
- `sparse::Bool`: if `false` call [`eigvals(A::QuantumObject; kwargs...)`](@ref), otherwise call [`eigsolve`](@ref). Default to `false`.
- `kwargs`: Additional keyword arguments passed to the solver. If `sparse=true`, the keyword arguments are passed to [`eigsolve`](@ref), otherwise to [`eigen`](@ref).

# Returns
- `::Vector{<:Number}`: a list of eigenvalues
"""
function eigenenergies(
    A::QuantumObject{OpType};
    sparse::Bool = false,
    kwargs...,
) where {OpType<:Union{Operator,SuperOperator}}
    if !sparse
        return eigvals(A; kwargs...)
    else
        return eigsolve(A; kwargs...).values
    end
end

@doc raw"""
    eigenstates(A::QuantumObject; sparse::Bool=false, kwargs...)

Calculate the eigenvalues and corresponding eigenvectors

# Arguments
- `A::QuantumObject`: the [`QuantumObject`](@ref) to solve eigenvalues and eigenvectors
- `sparse::Bool`: if `false` call [`eigen(A::QuantumObject; kwargs...)`](@ref), otherwise call [`eigsolve`](@ref). Default to `false`.
- `kwargs`: Additional keyword arguments passed to the solver. If `sparse=true`, the keyword arguments are passed to [`eigsolve`](@ref), otherwise to [`eigen`](@ref).

# Returns
- `::EigsolveResult`: containing the eigenvalues, the eigenvectors, and some information from the solver. see also [`EigsolveResult`](@ref)
"""
function eigenstates(
    A::QuantumObject{OpType};
    sparse::Bool = false,
    kwargs...,
) where {OpType<:Union{Operator,SuperOperator}}
    if !sparse
        return eigen(A; kwargs...)
    else
        return eigsolve(A; kwargs...)
    end
end
