using LinearAlgebra.BLAS: @blasfunc, BlasReal, BlasInt, BlasFloat, BlasComplex
using LinearAlgebra.BLAS: libblastrampoline
using LinearAlgebra: chkstride1, checksquare
using LinearAlgebra.LAPACK: chklapackerror
using Base: require_one_based_indexing
using ExponentialUtilities: arnoldi_step!

for (hseqr, elty) in
    ((:zhseqr_,:ComplexF64),
     (:chseqr_,:ComplexF32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOB, COMPZ
        #       INTEGER            N, ILO, IHI, LWORK, LDH, LDZ, INFO
        # *     ..
        # *     .. Array Arguments ..
        #       COMPLEX*16         H( LDH, * ), Z( LDZ, * ), WORK( * )
        function _hseqr!(job::AbstractChar, compz::AbstractChar, ilo::Integer, ihi::Integer,
                        H::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(H, Z)
            chkstride1(H)
            n = checksquare(H)
            checksquare(Z) == n || throw(DimensionMismatch())
            ldh = max(1, stride(H, 2))
            ldz = max(1, stride(Z, 2))
            w = similar(H, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hseqr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                    Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{BlasInt}),
                    job, compz, n, ilo, ihi,
                    H, ldh, w, Z, ldz, work,
                    lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            H, Z, w
        end
    end
end

for (hseqr, elty) in
    ((:dhseqr_,:Float64),
     (:shseqr_,:Float32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOB, COMPZ
        #       INTEGER            N, ILO, IHI, LWORK, LDH, LDZ, INFO
        # *     ..
        # *     .. Array Arguments ..
        #       COMPLEX*16         H( LDH, * ), Z( LDZ, * ), WORK( * )
        function _hseqr!(job::AbstractChar, compz::AbstractChar, ilo::Integer, ihi::Integer,
                        H::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(H, Z)
            chkstride1(H)
            n = checksquare(H)
            checksquare(Z) == n || throw(DimensionMismatch())
            ldh = max(1, stride(H, 2))
            ldz = max(1, stride(Z, 2))
            wr = similar(H, $elty, n)
            wi = similar(H, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hseqr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                    Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{BlasInt}),
                    job, compz, n, ilo, ihi,
                    H, ldh, wr, wi, Z, ldz, work,
                    lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            H, Z, complex.(wr, wi)
        end
    end
end
_hseqr!(H::StridedMatrix{T}, Z::StridedMatrix{T}) where {T<:BlasFloat} = _hseqr!('S', 'V', 1, size(H, 1), H, Z)
_hseqr!(H::StridedMatrix{T}) where {T<:BlasFloat} = _hseqr!('S', 'I', 1, size(H, 1), H, similar(H))

function _map_ldiv(linsolve, y, x)
    linsolve.b .= x
    y .= LinearSolve.solve!(linsolve).u
end

function _permuteschur!(T::AbstractMatrix{S}, Q::AbstractMatrix{S}, order::AbstractVector{<:Integer}) where {S<:BlasFloat}
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

function _eigsolve(A, b::AbstractVector, k::Integer = min(4, size(A, 1)), 
    m::Integer = min(10, size(A, 1)); tol::Real = 1e-8, maxiter::Integer = 200)

    Ks = ExponentialUtilities.arnoldi(A, b, m = m)
    V = Ks.V
    H = Ks.H
    f = ones(eltype(A), m)

    Vₘ = view(V, :, 1:m)
    Hₘ = view(H, 1:m, :)
    qₘ = view(V, :, m+1)
    βeₘ = view(H, m+1, :)
    β = H[m+1, m]
    Uₘ = one(Hₘ)
    
    Uₘᵥ = view(Uₘ, m, 1:m)

    cache0 = similar(b, m, m)
    cache1 = similar(b, size(V, 1), m)
    cache2 = similar(H, m)
    sorted_vals = Array{Int16}(undef, m)

    V₁ₖ = view(V, :, 1:k)
    Vₖ₊₁ = view(V, :, k+1)
    Hₖ₊₁₁ₖ = view(H, k+1, 1:k)
    cache1₁ₖ = view(cache1, :, 1:k)
    cache2₁ₖ = view(cache2, 1:k)

    M = typeof(cache0)

    numops = m
    iter = 1
    while iter < maxiter && count(x -> abs(x) < tol, f) < k
        # println( A * Vₘ ≈ Vₘ * M(Hₘ) + qₘ * M(transpose(βeₘ)) )     # SHOULD BE TRUE

        F = hessenberg!(Hₘ)
        copyto!(Uₘ, F.H.data)
        LAPACK.orghr!(1, m, Uₘ, F.τ)
        Tₘ, Uₘ, values = _hseqr!(F.H.data, Uₘ)
        sortperm!(sorted_vals, values, by = abs, rev = true)
        _permuteschur!(Tₘ, Uₘ, sorted_vals)
        mul!(f, Uₘᵥ, β)

        # println( A * Vₘ * Uₘ ≈ Vₘ * Uₘ * M(Tₘ) + qₘ * M(transpose(βeₘ)) * Uₘ )     # SHOULD BE TRUE

        copyto!(cache0, Uₘ)
        mul!(cache1, Vₘ, cache0)
        copyto!(V₁ₖ, cache1₁ₖ)
        copyto!(Vₖ₊₁, qₘ)
        mul!(cache2, transpose(Uₘ), βeₘ) # transpose(βeₘ) * Uₘ
        copyto!(Hₖ₊₁₁ₖ, cache2₁ₖ)

        # println( A * view(V, :, 1:k) ≈ view(V, :, 1:k) * M(view(H, 1:k, 1:k)) + qₘ * M(transpose(view(transpose(βeₘ) * Uₘ, 1:k))) )     # SHOULD BE TRUE

        for j in k+1:m
            β = arnoldi_step!(j, m, A, V, H, size(V, 1), 0)
        end

        numops += m-k-1
        iter+=1
    end

    F = hessenberg!(Hₘ)
    copyto!(Uₘ, F.H.data)
    LAPACK.orghr!(1, m, Uₘ, F.τ)
    Tₘ, Uₘ, values = _hseqr!(F.H.data, Uₘ)
    sortperm!(sorted_vals, values, by = abs, rev = true)
    _permuteschur!(Tₘ, Uₘ, sorted_vals)
    mul!(f, Uₘᵥ, β)

    vals = diag(view(Hₘ, 1:k, 1:k))
    select = Vector{BlasInt}(undef, 0)
    VR = LAPACK.trevc!('R', 'A', select, Tₘ)
    @inbounds for i in 1:size(VR, 2)
        normalize!(view(VR, :, i))
    end
    vecs = (Vₘ * M(Uₘ * VR))[:, 1:k]

    return vals, vecs, (iter, numops)
end

"""
    function eigsolve(A::QuantumObject; v0::Union{Nothing,AbstractVector}=nothing, 
        sigma::Union{Nothing, Real}=nothing, k::Integer = min(4, size(A, 1)), 
        krylovdim::Integer = min(10, size(A, 1)), tol::Real = 1e-8, maxiter::Integer = 200,
        solver::Union{Nothing, LinearSolve.SciMLLinearSolveAlgorithm} = nothing, showinfo::Bool=false, kwargs...)

Solve for the eigenvalues and eigenvectors of a matrix `A` using the Arnoldi method.
The keyword arguments are passed to the linear solver.
"""
function eigsolve(A::QuantumObject{<:AbstractMatrix}; v0::Union{Nothing,AbstractVector}=nothing, 
    sigma::Union{Nothing, Real}=nothing, k::Integer = min(4, size(A, 1)), 
    krylovdim::Integer = min(10, size(A, 1)), tol::Real = 1e-8, maxiter::Integer = 200,
    solver::Union{Nothing, LinearSolve.SciMLLinearSolveAlgorithm} = nothing, showinfo::Bool=false, kwargs...)

    return eigsolve(A.data; v0=v0, sigma=sigma, k=k, krylovdim=krylovdim, tol=tol, 
                    maxiter=maxiter, solver=solver, showinfo=showinfo, kwargs...)
end


function eigsolve(A::AbstractMatrix; v0::Union{Nothing,AbstractVector}=nothing, 
    sigma::Union{Nothing, Real}=nothing, k::Integer = min(4, size(A, 1)), 
    krylovdim::Integer = min(10, size(A, 1)), tol::Real = 1e-8, maxiter::Integer = 200,
    solver::Union{Nothing, LinearSolve.SciMLLinearSolveAlgorithm} = nothing, showinfo::Bool=false, kwargs...)

    T = eltype(A)
    isH = ishermitian(A)
    v0 === nothing && (v0 = normalize!(rand(eltype(T), size(A, 1))))

    if sigma === nothing
        vals, vecs, info = _eigsolve(A, v0, k, krylovdim, tol = tol, maxiter = maxiter)
        showinfo && println(info[1], " iterations, ", info[2], " applications of A")
    else
        Aₛ = A - sigma * I
        solver === nothing && (solver = isH ? KrylovJL_CG() : KrylovJL_GMRES())
        !haskey(kwargs, :Pl) && (kwargs = merge(kwargs, Dict(:Pl => ilu(Aₛ, τ=0.001))))

        prob = LinearProblem(Aₛ, v0)
        linsolve = init(prob, solver; kwargs...)
        Amap = LinearMap{T}((y,x) -> _map_ldiv(linsolve, y, x), length(v0))

        vals, vecs, info = _eigsolve(Amap, v0, k, krylovdim, tol = tol, maxiter = maxiter)
        showinfo && println(info[1], " iterations, ", info[2], " applications of A")
        @. vals = (1 + sigma * vals) / vals
    end

    isH && (vals = real.(vals))

    vals = vals[1:k]
    vecs = vecs[:, 1:k]

    return vals, vecs
end

"""
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

julia> E, U = eigen(H)
Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}
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

julia> ψ_1 = QuantumObject(U[:,1], dims=H.dims);

julia> expect(H, ψ_1) ≈ E[1]
true
```
"""
LinearAlgebra.eigen(A::QuantumObject{<:AbstractArray{T},OpType}; kwargs...) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = eigen(sparse_to_dense(A.data); kwargs...)

"""
    LinearAlgebra.eigvals(A::QuantumObject; kwargs...)

Same as [`eigen(A::QuantumObject; kwargs...)`](@ref) but for only the eigenvalues.
"""
LinearAlgebra.eigvals(A::QuantumObject{<:AbstractArray{T},OpType}; kwargs...) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = eigvals(sparse_to_dense(A.data); kwargs...)
