export WignerSolver, WignerClenshaw, WignerLaguerre, wigner

abstract type WignerSolver end

struct WignerClenshaw <: WignerSolver end

struct WignerLaguerre <: WignerSolver
    parallel::Bool
    tol::Float64
end

WignerLaguerre(; parallel = false, tol = 1.0e-14) = WignerLaguerre(parallel, tol)

@doc raw"""
    wigner(
        state::QuantumObject{OpType},
        xvec::AbstractVector,
        yvec::AbstractVector;
        g::Real = √2,
        method::WignerSolver = WignerClenshaw(),
    )

Generates the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) of `state` at points `xvec + 1im * yvec` in phase space. The `g` parameter is a scaling factor related to the value of ``\hbar`` in the commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2`` giving the default value ``\hbar=1``.

The `method` parameter can be either `WignerLaguerre()` or `WignerClenshaw()`. The former uses the Laguerre polynomial expansion of the Wigner function, while the latter uses the Clenshaw algorithm. The Laguerre expansion is faster for sparse matrices, while the Clenshaw algorithm is faster for dense matrices. The `WignerLaguerre` method has an optional `parallel` parameter which defaults to `true` and uses multithreading to speed up the calculation.

# Arguments
- `state::QuantumObject`: The quantum state for which the Wigner function is calculated. It can be either a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `xvec::AbstractVector`: The x-coordinates of the phase space grid.
- `yvec::AbstractVector`: The y-coordinates of the phase space grid.
- `g::Real`: The scaling factor related to the value of ``\hbar`` in the commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2``.
- `method::WignerSolver`: The method used to calculate the Wigner function. It can be either `WignerLaguerre()` or `WignerClenshaw()`, with `WignerClenshaw()` as default. The `WignerLaguerre` method has the optional `parallel` and `tol` parameters, with default values `true` and `1e-14`, respectively.

# Returns
- `W::Matrix`: The Wigner function of the state at the points `xvec + 1im * yvec` in phase space.

# Example
```jldoctest wigner
julia> ψ = fock(10, 0) + fock(10, 1) |> normalize

Quantum Object:   type=Ket()   dims=([10], [1])   size=(10,)
10-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im

julia> xvec = range(-5, 5, 200)
-5.0:0.05025125628140704:5.0

julia> wig = wigner(ψ, xvec, xvec);
```

or taking advantage of the parallel computation of the `WignerLaguerre` method

```jldoctest wigner
julia> ρ = ket2dm(ψ) |> to_sparse;

julia> wig = wigner(ρ, xvec, xvec, method=WignerLaguerre(parallel=true));

```
"""
function wigner(
        state::QuantumObject{OpType},
        xvec::AbstractVector,
        yvec::AbstractVector;
        g::Real = √2,
        method::WignerSolver = WignerClenshaw(),
    ) where {OpType <: Union{Bra, Ket, Operator}}
    ρ = ket2dm(state).data

    return _wigner(ρ, xvec, yvec, g, method)
end

function _wigner(
        ρ::AbstractArray,
        xvec::AbstractVector{T},
        yvec::AbstractVector{T},
        g::Real,
        method::WignerLaguerre,
    ) where {T <: Number}
    g = convert(T, g)
    X, Y = meshgrid(xvec, yvec)
    A = g / 2 * (X + 1im * Y)
    W = similar(A, T)
    W .= 0

    return _wigner_laguerre(ρ, A, W, g, method)
end

function _wigner(
        ρ::AbstractArray{T1},
        xvec::AbstractVector{T},
        yvec::AbstractVector{T},
        g::Real,
        method::WignerClenshaw,
    ) where {T1 <: Number, T <: Number}
    g = convert(T, g)
    M = size(ρ, 1)
    X, Y = meshgrid(xvec, yvec)
    A = g * (X + 1im * Y)

    B = abs.(A)
    B .*= B
    W = similar(A)
    W .= 2 * ρ[1, end]
    L = M - 1

    y0 = similar(B, T1)
    y1 = similar(B, T1)
    y0_old = copy(y0)
    res = similar(y0)

    while L > 0
        L -= 1
        ρdiag = _wig_laguerre_clenshaw!(res, L, B, lmul!(1 + Int(L != 0), diag(ρ, L)), y0, y1, y0_old)
        @. W = ρdiag + W * A / √(L + 1)
    end

    return @. real(W) * exp(-B / 2) * g^2 / 2 / π
end

function _wigner_laguerre(ρ::AbstractSparseArray, A::AbstractArray, W::AbstractArray, g::Real, method::WignerLaguerre)
    rows, cols, vals = findnz(ρ)
    B = @. 4 * abs2(A)

    if method.parallel
        iter = filter(x -> x[2] >= x[1], collect(zip(rows, cols, vals)))
        Wtot = similar(B, size(B)..., length(iter))
        Threads.@threads for i in eachindex(iter)
            m, n, ρmn = iter[i]
            m, n = m - 1, n - 1
            # Γ_mn = (1 + Int(m!=n)) * sqrt(gamma(m+1) / gamma(n+1))
            Γ_mn = (1 + Int(m != n)) * sqrt(exp(loggamma(m + 1) - loggamma(n + 1))) # Is this a good trick?
            Γ_mn = check_inf(Γ_mn)

            @. Wtot[:, :, i] = real(ρmn * (-1)^m * (2 * A)^(n - m) * Γ_mn * _genlaguerre(m, n - m, B))
        end
        W .= dropdims(sum(Wtot, dims = 3), dims = 3)
    else
        for i in Iterators.filter(x -> x[2] >= x[1], zip(rows, cols, vals))
            m, n, ρmn = i
            m, n = m - 1, n - 1
            # Γ_mn = (1 + Int(m!=n)) * sqrt(gamma(m+1) / gamma(n+1))
            Γ_mn = (1 + Int(m != n)) * sqrt(exp(loggamma(m + 1) - loggamma(n + 1))) # Is this a good trick?
            Γ_mn = check_inf(Γ_mn)

            @. W += real(ρmn * (-1)^m * (2 * A)^(n - m) * Γ_mn * _genlaguerre(m, n - m, B))
        end
    end

    return @. W * g^2 * exp(-B / 2) / 2 / π
end

function _wigner_laguerre(ρ::AbstractArray, A::AbstractArray, W::AbstractArray, g::Real, method::WignerLaguerre)
    tol = method.tol
    M = size(ρ, 1)
    B = @. 4 * abs2(A)

    if method.parallel
        throw(ArgumentError("Parallel version is not implemented for dense matrices"))
    else
        for m in 0:(M - 1)
            ρmn = ρ[m + 1, m + 1]
            abs(ρmn) > tol && (@. W += real(ρmn * (-1)^m * _genlaguerre(m, 0, B)))
            for n in (m + 1):(M - 1)
                ρmn = ρ[m + 1, n + 1]
                # Γ_mn = sqrt(gamma(m+1) / gamma(n+1))
                Γ_mn = sqrt(exp(loggamma(m + 1) - loggamma(n + 1))) # Is this a good trick?
                Γ_mn = check_inf(Γ_mn)

                abs(ρmn) > tol && (@. W += 2 * real(ρmn * (-1)^m * (2 * A)^(n - m) * Γ_mn * _genlaguerre(m, n - m, B)))
            end
        end
    end

    return @. W * g^2 * exp(-B / 2) / 2 / π
end

# function _genlaguerre(n::Int, α::Int, x::T) where {T<:Number}
#     t = binomial(n+α,n)
#     L = t
#     for k = 1:n
#         t *= (-x) * (n-k+1) /  (k * (k+α))
#         L += t
#     end
#     return L
# end
# This is a little bit slower, but it supports GPU
function _genlaguerre(n::Int, α::Number, x::T) where {T <: Number}
    α = convert(T, α)
    p0, p1 = one(T), -x + (α + 1)
    n == 0 && return p0
    for k in 1:(n - 1)
        p1, p0 = ((2k + α + 1) / (k + 1) - x / (k + 1)) * p1 - (k + α) / (k + 1) * p0, p1
    end
    return p1
end

check_inf(x::T) where {T} =
if isinf(x)
    return x > zero(T) ? floatmax(T) : -floatmax(T)
else
    return x
end

function _wig_laguerre_clenshaw!(res, L::Int, x, c, y0, y1, y0_old)
    length(c) == 1 && return c[1]
    length(c) == 2 && return @. c[1] - c[2] * (L + 1 - x) / sqrt(L + 1)

    y0 .= c[end - 1]
    y1 .= c[end]

    k = length(c)
    for i in range(3, length(c), step = 1)
        k -= 1
        copyto!(y0_old, y0)
        @. y0 = c[end + 1 - i] - y1 * sqrt((k - 1) * (L + k - 1) / ((L + k) * k))
        @. y1 = y0_old - y1 * ((L + 2 * k - 1) - x) / sqrt((L + k) * k)
    end

    @. res = y0 - y1 * (L + 1 - x) / sqrt(L + 1)

    return res
end
