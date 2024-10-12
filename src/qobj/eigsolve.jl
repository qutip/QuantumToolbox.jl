#=
Eigen solvers and results for QuantumObject
=#

export EigsolveResult
export eigenenergies, eigenstates, eigsolve
export eigsolve_al

@doc raw"""
    struct EigsolveResult{T1<:Vector{<:Number}, T2<:AbstractMatrix{<:Number}, ObjType<:Union{Nothing,OperatorQuantumObject,SuperOperatorQuantumObject},N}
        values::T1
        vectors::T2
        type::ObjType
        dims::SVector{N,Int}
        iter::Int
        numops::Int
        converged::Bool
    end

A struct containing the eigenvalues, the eigenvectors, and some information from the solver

# Fields
- `values::AbstractVector`: the eigenvalues
- `vectors::AbstractMatrix`: the transformation matrix (eigenvectors)
- `type::Union{Nothing,QuantumObjectType}`: the type of [`QuantumObject`](@ref), or `nothing` means solving eigen equation for general matrix
- `dims::SVector`: the `dims` of [`QuantumObject`](@ref)
- `iter::Int`: the number of iteration during the solving process
- `numops::Int` : number of times the linear map was applied in krylov methods
- `converged::Bool`: Whether the result is converged

# Examples
One can obtain the eigenvalues and the corresponding [`QuantumObject`](@ref)-type eigenvectors by:
```
julia> result = eigenstates(sigmax());

julia> λ, ψ, T = result;

julia> λ
2-element Vector{ComplexF64}:
 -1.0 + 0.0im
  1.0 + 0.0im

julia> ψ
2-element Vector{QuantumObject{Vector{ComplexF64}, KetQuantumObject}}:
 QuantumObject{Vector{ComplexF64}, KetQuantumObject}(ComplexF64[-0.7071067811865475 + 0.0im, 0.7071067811865475 + 0.0im], KetQuantumObject(), [2])
 QuantumObject{Vector{ComplexF64}, KetQuantumObject}(ComplexF64[0.7071067811865475 + 0.0im, 0.7071067811865475 + 0.0im], KetQuantumObject(), [2])

julia> T
2×2 Matrix{ComplexF64}:
 -0.707107+0.0im  0.707107+0.0im
  0.707107+0.0im  0.707107+0.0im
```
"""
struct EigsolveResult{
    T1<:Vector{<:Number},
    T2<:AbstractMatrix{<:Number},
    ObjType<:Union{Nothing,OperatorQuantumObject,SuperOperatorQuantumObject},
    N,
}
    values::T1
    vectors::T2
    type::ObjType
    dims::SVector{N,Int}
    iter::Int
    numops::Int
    converged::Bool
end

Base.iterate(res::EigsolveResult) = (res.values, Val(:vector_list))
Base.iterate(res::EigsolveResult{T1,T2,Nothing}, ::Val{:vector_list}) where {T1,T2} =
    ([res.vectors[:, k] for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1,T2,OperatorQuantumObject}, ::Val{:vector_list}) where {T1,T2} =
    ([QuantumObject(res.vectors[:, k], Ket, res.dims) for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult{T1,T2,SuperOperatorQuantumObject}, ::Val{:vector_list}) where {T1,T2} =
    ([QuantumObject(res.vectors[:, k], OperatorKet, res.dims) for k in 1:length(res.values)], Val(:vectors))
Base.iterate(res::EigsolveResult, ::Val{:vectors}) = (res.vectors, Val(:done))
Base.iterate(res::EigsolveResult, ::Val{:done}) = nothing

function Base.show(io::IO, res::EigsolveResult)
    println(io, "EigsolveResult:   type=", res.type, "   dims=", res.dims)
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

function _permuteschur!(
    T::AbstractMatrix{S},
    Q::AbstractMatrix{S},
    order::AbstractVector{<:Integer},
) where {S<:BlasFloat}
    n = checksquare(T)
    p = collect(order) # makes copy cause will be overwritten
    @inbounds for i in eachindex(p)
        ifirst::BlasInt = p[i]
        ilast::BlasInt = i
        LAPACK.trexc!(ifirst, ilast, T, Q)
        for k in (i+1):length(p)
            if p[k] < p[i]
                p[k] += 1
            end
        end
    end
    return T, Q
end

function _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, m, β, sorted_vals)
    F = hessenberg!(Hₘ)
    copyto!(Uₘ, Hₘ)
    LAPACK.orghr!(1, m, Uₘ, F.τ)
    Tₘ, Uₘ, values = hseqr!(Hₘ, Uₘ)
    sortperm!(sorted_vals, values, by = abs, rev = true)
    _permuteschur!(Tₘ, Uₘ, sorted_vals)
    mul!(f, Uₘᵥ, β)

    return Tₘ, Uₘ
end

function _eigsolve(
    A,
    b::AbstractVector{T},
    type::ObjType,
    dims::SVector,
    k::Int = 1,
    m::Int = max(20, 2 * k + 1);
    tol::Real = 1e-8,
    maxiter::Int = 200,
) where {T<:BlasFloat,ObjType<:Union{Nothing,OperatorQuantumObject,SuperOperatorQuantumObject}}
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
    sorted_vals = Array{Int16}(undef, m)

    V₁ₖ = view(V, :, 1:k)
    Vₖ₊₁ = view(V, :, k + 1)
    Hₖ₊₁₁ₖ = view(H, k + 1, 1:k)
    cache1₁ₖ = view(cache1, :, 1:k)
    cache2₁ₖ = view(cache2, 1:k)

    M = typeof(cache0)

    Tₘ, Uₘ = _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, m, β, sorted_vals)

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

        for j in k+1:m
            β = arnoldi_step!(A, V, H, j)
            if β < tol
                numops += j - k - 1
                break
            end
        end

        # println( A * Vₘ ≈ Vₘ * M(Hₘ) + qₘ * M(transpose(βeₘ)) )     # SHOULD BE TRUE

        Tₘ, Uₘ = _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, m, β, sorted_vals)

        numops += m - k - 1
        iter += 1
    end

    vals = diag(view(Tₘ, 1:k, 1:k))
    select = Vector{BlasInt}(undef, 0)
    VR = LAPACK.trevc!('R', 'A', select, Tₘ)
    @inbounds for i in 1:size(VR, 2)
        normalize!(view(VR, :, i))
    end
    mul!(cache1, Vₘ, M(Uₘ * VR))
    vecs = cache1[:, 1:k]

    return EigsolveResult(vals, vecs, type, dims, iter, numops, (iter < maxiter))
end

@doc raw"""
    eigsolve(A::QuantumObject; 
        v0::Union{Nothing,AbstractVector}=nothing, 
        sigma::Union{Nothing, Real}=nothing,
        k::Int = 1,
        krylovdim::Int = max(20, 2*k+1),
        tol::Real = 1e-8,
        maxiter::Int = 200,
        solver::Union{Nothing, SciMLLinearSolveAlgorithm} = nothing,
        kwargs...)

Solve for the eigenvalues and eigenvectors of a matrix `A` using the Arnoldi method.

# Notes
- For more details about `solver` and extra `kwargs`, please refer to [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/)

# Returns
- `EigsolveResult`: A struct containing the eigenvalues, the eigenvectors, and some information about the eigsolver
"""
function eigsolve(
    A::QuantumObject{<:AbstractMatrix};
    v0::Union{Nothing,AbstractVector} = nothing,
    sigma::Union{Nothing,Real} = nothing,
    k::Int = 1,
    krylovdim::Int = max(20, 2 * k + 1),
    tol::Real = 1e-8,
    maxiter::Int = 200,
    solver::Union{Nothing,SciMLLinearSolveAlgorithm} = nothing,
    kwargs...,
)
    return eigsolve(
        A.data;
        v0 = v0,
        type = A.type,
        dims = A.dims,
        sigma = sigma,
        k = k,
        krylovdim = krylovdim,
        tol = tol,
        maxiter = maxiter,
        solver = solver,
        kwargs...,
    )
end

function eigsolve(
    A;
    v0::Union{Nothing,AbstractVector} = nothing,
    type::Union{Nothing,OperatorQuantumObject,SuperOperatorQuantumObject} = nothing,
    dims = SVector{0,Int}(),
    sigma::Union{Nothing,Real} = nothing,
    k::Int = 1,
    krylovdim::Int = max(20, 2 * k + 1),
    tol::Real = 1e-8,
    maxiter::Int = 200,
    solver::Union{Nothing,SciMLLinearSolveAlgorithm} = nothing,
    kwargs...,
)
    T = eltype(A)
    isH = ishermitian(A)
    v0 === nothing && (v0 = normalize!(rand(T, size(A, 1))))

    dims = SVector(dims)

    if sigma === nothing
        res = _eigsolve(A, v0, type, dims, k, krylovdim, tol = tol, maxiter = maxiter)
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

        res = _eigsolve(Amap, v0, type, dims, k, krylovdim, tol = tol, maxiter = maxiter)
        vals = @. (1 + sigma * res.values) / res.values
    end

    return EigsolveResult(vals, res.vectors, res.type, res.dims, res.iter, res.numops, res.converged)
end

@doc raw"""
    eigsolve_al(H::QuantumObject,
        T::Real, c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        H_t::Union{Nothing,Function}=nothing,
        params::NamedTuple=NamedTuple(),
        ρ0::Union{Nothing, AbstractMatrix} = nothing,
        k::Int=1,
        krylovdim::Int=min(10, size(H, 1)),
        maxiter::Int=200,
        eigstol::Real=1e-6,
        kwargs...)

Solve the eigenvalue problem for a Liouvillian superoperator `L` using the Arnoldi-Lindblad method.

# Arguments
- `H`: The Hamiltonian (or directly the Liouvillian) of the system.
- `T`: The time at which to evaluate the time evolution
- `c_ops`: A vector of collapse operators. Default is `nothing` meaning the system is closed.
- `alg`: The differential equation solver algorithm
- `H_t`: A function `H_t(t)` that returns the additional term at time `t`
- `params`: A dictionary of additional parameters
- `ρ0`: The initial density matrix. If not specified, a random density matrix is used
- `k`: The number of eigenvalues to compute
- `krylovdim`: The dimension of the Krylov subspace
- `maxiter`: The maximum number of iterations for the eigsolver
- `eigstol`: The tolerance for the eigsolver
- `kwargs`: Additional keyword arguments passed to the differential equation solver

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- `EigsolveResult`: A struct containing the eigenvalues, the eigenvectors, and some information about the eigsolver

# References
- [1] Minganti, F., & Huybrechts, D. (2022). Arnoldi-Lindblad time evolution: Faster-than-the-clock algorithm for the spectrum of time-independent and Floquet open quantum systems. Quantum, 6, 649.
"""
function eigsolve_al(
    H::QuantumObject{MT1,HOpType},
    T::Real,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    H_t::Union{Nothing,Function} = nothing,
    params::NamedTuple = NamedTuple(),
    ρ0::AbstractMatrix = rand_dm(prod(H.dims)).data,
    k::Int = 1,
    krylovdim::Int = min(10, size(H, 1)),
    maxiter::Int = 200,
    eigstol::Real = 1e-6,
    kwargs...,
) where {MT1<:AbstractMatrix,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    L = liouvillian(H, c_ops)
    prob = mesolveProblem(
        L,
        QuantumObject(ρ0, type=Operator, dims = H.dims),
        [zero(T), T];
        alg = alg,
        H_t = H_t,
        params = params,
        progress_bar = Val(false),
        kwargs...,
    )
    integrator = init(prob, alg)

    # prog = ProgressUnknown(desc="Applications:", showspeed = true, enabled=progress)

    Lmap = ArnoldiLindbladIntegratorMap(eltype(MT1), size(L), integrator)

    res = _eigsolve(Lmap, mat2vec(ρ0), L.type, L.dims, k, krylovdim, maxiter = maxiter, tol = eigstol)
    # finish!(prog)

    vals = similar(res.values)
    vecs = similar(res.vectors)

    for i in eachindex(res.values)
        vec = view(res.vectors, :, i)
        vals[i] = dot(vec, L.data, vec)
        @. vecs[:, i] = vec * exp(-1im * angle(vec[1]))
    end

    return EigsolveResult(vals, vecs, res.type, res.dims, res.iter, res.numops, res.converged)
end

@doc raw"""
    LinearAlgebra.eigen(A::QuantumObject; kwargs...)

Calculates the eigenvalues and eigenvectors of the [`QuantumObject`](@ref) `A` using
the Julia [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) package.

```
julia> a = destroy(5);

julia> H = a + a'
Quantum Object:   type=Operator   dims=[5]   size=(5, 5)   ishermitian=true
5×5 SparseMatrixCSC{ComplexF64, Int64} with 8 stored entries:
     ⋅          1.0+0.0im          ⋅              ⋅          ⋅
 1.0+0.0im          ⋅      1.41421+0.0im          ⋅          ⋅
     ⋅      1.41421+0.0im          ⋅      1.73205+0.0im      ⋅
     ⋅              ⋅      1.73205+0.0im          ⋅      2.0+0.0im
     ⋅              ⋅              ⋅          2.0+0.0im      ⋅

julia> E, ψ, U = eigen(H)
EigsolveResult:   type=Operator   dims=[5]
values:
5-element Vector{Float64}:
 -2.8569700138728
 -1.3556261799742608
  1.3322676295501878e-15
  1.3556261799742677
  2.8569700138728056
vectors:
5×5 Matrix{ComplexF64}:
  0.106101+0.0im  -0.471249-0.0im  …   0.471249-0.0im  0.106101-0.0im
 -0.303127-0.0im   0.638838+0.0im      0.638838+0.0im  0.303127-0.0im
  0.537348+0.0im  -0.279149-0.0im      0.279149-0.0im  0.537348-0.0im
 -0.638838-0.0im  -0.303127-0.0im     -0.303127-0.0im  0.638838+0.0im
  0.447214+0.0im   0.447214+0.0im     -0.447214-0.0im  0.447214-0.0im

julia> expect(H, ψ[1]) ≈ E[1]
true
```
"""
function LinearAlgebra.eigen(
    A::QuantumObject{MT,OpType};
    kwargs...,
) where {MT<:AbstractMatrix,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    F = eigen(sparse_to_dense(A.data); kwargs...)
    # This fixes a type inference issue. But doesn't work for GPU arrays
    E::mat2vec(sparse_to_dense(MT)) = F.values
    U::sparse_to_dense(MT) = F.vectors

    return EigsolveResult(E, U, A.type, A.dims, 0, 0, true)
end

@doc raw"""
    LinearAlgebra.eigvals(A::QuantumObject; kwargs...)

Same as [`eigen(A::QuantumObject; kwargs...)`](@ref) but for only the eigenvalues.
"""
LinearAlgebra.eigvals(
    A::QuantumObject{<:AbstractArray{T},OpType};
    kwargs...,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    eigvals(sparse_to_dense(A.data); kwargs...)

@doc raw"""
    eigenenergies(A::QuantumObject; sparse::Bool=false, kwargs...)

Calculate the eigenenergies

# Arguments
- `A::QuantumObject`: the [`QuantumObject`](@ref) to solve eigenvalues
- `sparse::Bool`: if `false` call [`eigvals(A::QuantumObject; kwargs...)`](@ref), otherwise call [`eigsolve`](@ref). Default to `false`.
- `kwargs`: Additional keyword arguments passed to the solver

# Returns
- `::Vector{<:Number}`: a list of eigenvalues
"""
function eigenenergies(
    A::QuantumObject{<:AbstractArray{T},OpType};
    sparse::Bool = false,
    kwargs...,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
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
- `kwargs`: Additional keyword arguments passed to the solver

# Returns
- `::EigsolveResult`: containing the eigenvalues, the eigenvectors, and some information from the solver. see also [`EigsolveResult`](@ref)
"""
function eigenstates(
    A::QuantumObject{<:AbstractArray{T},OpType};
    sparse::Bool = false,
    kwargs...,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    if !sparse
        return eigen(A; kwargs...)
    else
        return eigsolve(A; kwargs...)
    end
end
