abstract type WignerSolver end
struct WignerClenshaw <: WignerSolver end
struct WignerLaguerre <: WignerSolver
    parallel::Bool
    tol::Float64
end

WignerLaguerre(;parallel=false, tol=1e-14) = WignerLaguerre(parallel, tol)

@doc raw"""
    wigner(state::QuantumObject, xvec::AbstractVector, yvec::AbstractVector; g::Real=√2,
        solver::WignerSolver=WignerLaguerre())

Generates the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution)
of `state` at points `xvec + 1im * yvec`. The `g` parameter is a scaling factor related to the value of ``\hbar`` in the
commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2`` giving the default value ``\hbar=1``.

The `solver` parameter can be either `WignerLaguerre()` or `WignerClenshaw()`. The former uses the Laguerre polynomial
expansion of the Wigner function, while the latter uses the Clenshaw algorithm. The Laguerre expansion is faster for
sparse matrices, while the Clenshaw algorithm is faster for dense matrices. The `WignerLaguerre` solver has an optional
`parallel` parameter which defaults to `true` and uses multithreading to speed up the calculation.
"""
function wigner(state::QuantumObject{<:AbstractArray{T},OpType}, xvec::AbstractVector,
    yvec::AbstractVector; g::Real=√2, solver::MySolver=WignerLaguerre()) where 
    {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject},
    MySolver<:WignerSolver}

    if isket(state)
        ρ = (state * state').data
    elseif isbra(state)
        ρ = (state' * state).data
    else
        ρ = state.data
    end

    return _wigner(ρ, xvec, yvec, g, solver)
end

function _wigner(ρ::AbstractArray, xvec::AbstractVector, yvec::AbstractVector,
    g::Real, solver::WignerLaguerre)
    
    X, Y = meshgrid(xvec, yvec)
    A = g / 2 * (X + 1im * Y)
    W = similar(A, Float64)
    W .= 0

    return _wigner_laguerre(ρ, A, W, g, solver)
end

function _wigner(ρ::AbstractArray, xvec::AbstractVector, yvec::AbstractVector,
    g::Real, solver::WignerClenshaw)
    
    M = size(ρ, 1)
    X, Y = meshgrid(xvec, yvec)
    A = g * (X + 1im * Y)

    B = abs.(A)
    B .*= B
    W = similar(A)
    W .= 2 * ρ[1, end]
    L = M - 1

    while L > 0
        L -= 1
        ρdiag = (L == 0) ? _wig_laguerre_clenshaw(L, B, diag(ρ, L)) : _wig_laguerre_clenshaw(L, B, 2*diag(ρ, L))
        @. W = ρdiag + W * A / √(L + 1)
    end

    return @. real(W) * exp(-B / 2) * (g^2 / (2π))
end

function _wigner_laguerre(ρ::AbstractSparseArray, A::AbstractArray, W::AbstractArray, g::Real, solver::WignerLaguerre)
    rows, cols, vals = findnz(ρ)
    B = @. 4 * abs2(A)

    if solver.parallel
        iter = filter(x->x[2]>=x[1], collect(zip(rows, cols, vals)))
        Wtot = similar(B, size(B)..., length(iter))
        Threads.@threads for i in eachindex(iter)
            m, n, ρmn = iter[i]
            m, n = m-1, n-1

            if m==n
                @. Wtot[:,:,i] = real(ρmn * (-1)^m * _genlaguerre(m, 0, B))
            else
                @. Wtot[:,:,i] = 2 * real(ρmn * (-1)^m * (2 * A)^(n - m) * sqrt(factorial(big(m)) / factorial(big(n))) *
                     _genlaguerre(m, n - m, B))
            end
        end
        W .= dropdims(sum(Wtot, dims=3), dims=3)
    else
        for i in Iterators.filter(x->x[2]>=x[1], zip(rows, cols, vals))
            m, n, ρmn = i
            m, n = m-1, n-1

            if m==n
                @. W += real(ρmn * (-1)^m * _genlaguerre(m, 0, B))
            else
                @. W += 2 * real(ρmn * (-1)^m * (2 * A)^(n - m) * sqrt(factorial(big(m)) / factorial(big(n))) *
                     _genlaguerre(m, n - m, B))
            end
        end
    end

    return @. W * g^2 * exp(-B / 2) / (2π)
end

function _wigner_laguerre(ρ::AbstractArray, A::AbstractArray, W::AbstractArray, g::Real, solver::WignerLaguerre)
    tol = solver.tol
    M = size(ρ, 1)
    B = @. 4 * abs2(A)
    
    if solver.parallel
        throw(ArgumentError("Parallel version is not implemented for dense matrices"))
    else
        for m in 0:M-1
            ρmn = ρ[m+1, m+1]
            abs(ρmn) > tol && (@. W += real(ρmn * (-1)^m * _genlaguerre(m, 0, B)))
            for n in m+1:M-1
                ρmn = ρ[m+1, n+1]
                abs(ρmn) > tol && (@. W += 2 * real(ρmn * (-1)^m * (2 * A)^(n - m) * sqrt(factorial(big(m)) / factorial(big(n))) *
                     _genlaguerre(m, n - m, B)))
            end
        end
    end

    return @. W * g^2 * exp(-B / 2) / (2π)
end

_genlaguerre(n::Integer, α::Integer, x::T) where {T} = binomial(n+α,n) * HypergeometricFunctions.M(-n, α+1, x)

function _wig_laguerre_clenshaw(L, x, c)
    if length(c) == 1
        y0 = c[1]
        y1 = 0
    elseif length(c) == 2
        y0 = c[1]
        y1 = c[2]
    else
        k = length(c)
        y0 = c[end-1]
        y1 = c[end]
        for i in range(3, length(c), step=1)
            k -= 1
            y0, y1 = @. c[end+1-i] - y1 * ((k - 1) * (L + k - 1) / ((L + k) * k))^0.5, y0 - y1 * ((L + 2 * k - 1) - x) * ((L + k) * k)^(-0.5)
        end
    end
    return @. y0 - y1 * ((L + 1) - x) * (L + 1)^(-0.5)
end