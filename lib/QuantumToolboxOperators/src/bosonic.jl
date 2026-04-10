# ──────────────────────────────────────────────────────────────────────────────
# Bosonic matrix-free quantum operators
# ──────────────────────────────────────────────────────────────────────────────

"""
    BosonicOperator{T} <: AbstractSciMLOperator{T}

Abstract supertype for all matrix-free bosonic operators.
"""
abstract type BosonicOperator{T} <: AbstractSciMLOperator{T} end

# ─── Type-aware coefficient helpers ──────────────────────────────────────────
# Avoid sqrt(::Int) → Float64 promotion when working with Float32 / other types.

_sqrt_coeff(::Type{T}, i::Int) where {T} = sqrt(real(T)(i))

function _power_coeff(::Type{T}, i::Int, k::Int) where {T}
    c = one(T)
    for m in i:(i + k - 1)
        c *= real(T)(m)
    end
    return sqrt(c)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyOperator{T,Precomp} <: AbstractSciMLOperator{T}

Matrix-free bosonic annihilation (destroy) operator.

Action on a Fock-basis vector: ``â |n⟩ = √n |n-1⟩``.

The creation operator ``â†`` is obtained via `adjoint(a)`, which returns
an `AdjointOperator{T, DestroyOperator}` — no separate type is needed.

When `Precomp = true`, the coefficients ``√1, √2, …, √(N-1)`` are precomputed
and stored. When `Precomp = false`, they are computed on the fly.
"""
struct DestroyOperator{T, Precomp, VT} <: BosonicOperator{T}
    N::Int
    coeffs::VT

    function DestroyOperator{T, true}(N::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        coeffs = [_sqrt_coeff(T, i) for i in 1:(N - 1)]
        return new{T, true, typeof(coeffs)}(N, coeffs)
    end

    function DestroyOperator{T, false}(N::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T, false, Nothing}(N, nothing)
    end

    function DestroyOperator(N::Int, coeffs::Union{Nothing, AbstractVector})
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        if isnothing(coeffs)
            return new{Float64, false, Nothing}(N, nothing)
        else
            return new{eltype(coeffs), true, typeof(coeffs)}(N, coeffs)
        end
    end
end

DestroyOperator{T}(N::Int) where {T} = DestroyOperator{T, true}(N)
DestroyOperator(N::Int) = DestroyOperator{Float64}(N)

Base.size(L::DestroyOperator) = (L.N, L.N)
Base.size(L::DestroyOperator, n::Int) = L.N

islinear(::DestroyOperator) = true
has_adjoint(::DestroyOperator) = true

# ─── mul! with precomputed coefficients ──────────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyOperator{T, true}, v::AbstractVecOrMat) where {T}
    N = L.N

    # @views w[N:N] .= zero(eltype(w)) # This doesn't work on Reactant
    fill!(w, zero(eltype(w)))
    @views w[1:(N - 1)] .= L.coeffs .* v[2:N]
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyOperator{T, true}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.N
    # @views w[N:N] .= β .* w[N:N] # This doesn't work on Reactant
    lmul!(β, w)
    @views w[1:(N - 1)] .= α .* L.coeffs .* v[2:N] .+ β .* w[1:(N - 1)]
    return w
end

# ─── mul! without precomputed coefficients ───────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyOperator{T, false}, v::AbstractVecOrMat) where {T}
    N = L.N
    # @inbounds @simd for i in 1:(N - 1)
    #     w[i] = _sqrt_coeff(T, i) * v[i + 1]
    # end
    # @inbounds w[N] = zero(eltype(w))

    fill!(w, zero(eltype(w)))
    @views w[1:(N - 1)] .= _sqrt_coeff.(T, 1:(N - 1)) .* v[2:N]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyOperator{T, false}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.N
    # @inbounds @simd for i in 1:(N - 1)
    #     w[i] = α * _sqrt_coeff(T, i) * v[i + 1] + β * w[i]
    # end
    # @inbounds w[N] = β * w[N]

    lmul!(β, w)
    @views w[1:(N - 1)] .= α .* _sqrt_coeff.(T, 1:(N - 1)) .* v[2:N] .+ β .* w[1:(N - 1)]

    return w
end

# ─── Adjoint: creation operator â† ───────────────────────────────────────────
# â† |n⟩ = √(n+1) |n+1⟩  ⟹  w[1] = 0;  w[i+1] = √i · v[i]  for i = 1…N-1

const AdjointDestroyOperator{T, Precomp} = AdjointOperator{T, <:DestroyOperator{T, Precomp}} where {T, Precomp}

# ─── Adjoint mul! with precomputed coefficients ─────────────────────────────

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T, true}}, v::AbstractVecOrMat,
    ) where {T}
    N = L.L.N
    # @views w[1:1] .= zero(eltype(w)) # This doesn't work on Reactant
    fill!(w, zero(eltype(w)))
    @views w[2:N] .= L.L.coeffs .* v[1:(N - 1)]
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T, true}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.L.N
    # @views w[1:1] .= β .* w[1:1] # This doesn't work on Reactant
    lmul!(β, w)
    @views w[2:N] .= α .* L.L.coeffs .* v[1:(N - 1)] .+ β .* w[2:N]
    return w
end

# ─── Adjoint mul! without precomputed coefficients ──────────────────────────

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T, false}}, v::AbstractVecOrMat,
    ) where {T}
    N = L.L.N

    # @inbounds w[1] = zero(eltype(w))
    # @inbounds @simd for i in 1:(N - 1)
    #     w[i + 1] = _sqrt_coeff(T, i) * v[i]
    # end

    fill!(w, zero(eltype(w)))
    @views w[2:N] .= _sqrt_coeff.(T, 1:(N - 1)) .* v[1:(N - 1)]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T, false}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.L.N

    # @inbounds w[1] = β * w[1]
    # @inbounds @simd for i in 1:(N - 1)
    #     w[i + 1] = α * _sqrt_coeff(T, i) * v[i] + β * w[i + 1]
    # end

    lmul!(β, w)
    @views w[2:N] .= α .* _sqrt_coeff.(T, 1:(N - 1)) .* v[1:(N - 1)] .+ β .* w[2:N]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NumberOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    NumberOperator{T,Precomp} <: AbstractSciMLOperator{T}

Matrix-free bosonic number operator with optional shift.

Action: ``w_i = (i - 1 + \\text{shift}) \\cdot v_i``.

  - `shift = 0` represents ``\\hat{n} = â†â`` with eigenvalues ``0, 1, 2, …``
  - `shift = 1` represents ``ââ†`` with eigenvalues ``1, 2, 3, …``

When `Precomp = true`, the diagonal coefficients are precomputed and stored.
When `Precomp = false`, they are computed on the fly.
"""
struct NumberOperator{T, Precomp, VT} <: BosonicOperator{T}
    N::Int
    shift::Int
    coeffs::VT

    function NumberOperator{T, true}(N::Int; shift::Int = 0) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        coeffs = [real(T)(i - 1 + shift) for i in 1:N]
        return new{T, true, typeof(coeffs)}(N, shift, coeffs)
    end

    function NumberOperator{T, false}(N::Int; shift::Int = 0) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T, false, Nothing}(N, shift, nothing)
    end

    function NumberOperator(N::Int, shift::Int, coeffs::Union{Nothing, AbstractVector})
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        if isnothing(coeffs)
            return new{Float64, false, Nothing}(N, shift, nothing)
        else
            return new{eltype(coeffs), true, typeof(coeffs)}(N, shift, coeffs)
        end
    end
end

NumberOperator{T}(N::Int; shift::Int = 0) where {T} = NumberOperator{T, true}(N; shift)
NumberOperator(N::Int; shift::Int = 0) = NumberOperator{Float64}(N; shift)

Base.size(L::NumberOperator) = (L.N, L.N)
Base.size(L::NumberOperator, n::Int) = L.N

islinear(::NumberOperator) = true
has_adjoint(::NumberOperator) = true

# NumberOperator is Hermitian: adjoint returns itself
Base.adjoint(L::NumberOperator) = L

# ─── mul! with precomputed coefficients ──────────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NumberOperator{T, true}, v::AbstractVecOrMat) where {T}
    w .= L.coeffs .* v
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NumberOperator{T, true}, v::AbstractVecOrMat, α, β,
    ) where {T}
    w .= α .* L.coeffs .* v .+ β .* w
    return w
end

# ─── mul! without precomputed coefficients ───────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NumberOperator{T, false}, v::AbstractVecOrMat) where {T}
    N, shift = L.N, L.shift

    # @inbounds @simd for i in 1:N
    #     w[i] = T(i - 1 + shift) * v[i]
    # end

    fill!(w, zero(eltype(w)))
    w .= real(T).(shift:shift+N-1) .* v

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NumberOperator{T, false}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, shift = L.N, L.shift

    # @inbounds @simd for i in 1:N
    #     w[i] = α * T(i - 1 + shift) * v[i] + β * w[i]
    # end

    lmul!(β, w)
    w .= α .* real(T).(shift:shift+N-1) .* v .+ β .* w

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyPowerOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyPowerOperator{T,Precomp} <: AbstractSciMLOperator{T}

Matrix-free ``â^k`` (k-th power of the bosonic annihilation operator).

Action on Fock basis: ``â^k |n⟩ = \\sqrt{n! / (n-k)!}\\, |n-k⟩`` for ``n ≥ k``, else ``0``.

When `Precomp = true`, the coefficients are precomputed and stored.
When `Precomp = false`, they are computed on the fly.
"""
struct DestroyPowerOperator{T, Precomp, VT} <: BosonicOperator{T}
    N::Int
    k::Int
    coeffs::VT

    function DestroyPowerOperator{T, true}(N::Int, k::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        coeffs = [_power_coeff(T, i, k) for i in 1:(N - k)]
        return new{T, true, typeof(coeffs)}(N, k, coeffs)
    end

    function DestroyPowerOperator{T, false}(N::Int, k::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        return new{T, false, Nothing}(N, k, nothing)
    end

    function DestroyPowerOperator(N::Int, k::Int, coeffs::Union{Nothing, AbstractVector})
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        if isnothing(coeffs)
            return new{Float64, false, Nothing}(N, k, nothing)
        else
            return new{eltype(coeffs), true, typeof(coeffs)}(N, k, coeffs)
        end
    end
end

DestroyPowerOperator{T}(N::Int, k::Int) where {T} = DestroyPowerOperator{T, true}(N, k)
DestroyPowerOperator(N::Int, k::Int) = DestroyPowerOperator{Float64}(N, k)

Base.size(L::DestroyPowerOperator) = (L.N, L.N)
Base.size(L::DestroyPowerOperator, n::Int) = L.N

islinear(::DestroyPowerOperator) = true
has_adjoint(::DestroyPowerOperator) = true

# ─── mul! with precomputed coefficients ──────────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyPowerOperator{T, true}, v::AbstractVecOrMat) where {T}
    N, k = L.N, L.k

    # @views w[(N - k + 1):N] .= zero(eltype(w)) # This doesn't work on Reactant
    fill!(w, zero(eltype(w)))
    @views w[1:(N - k)] .= L.coeffs .* v[(k + 1):N]
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyPowerOperator{T, true}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.N, L.k

    # @views w[(N - k + 1):N] .= β .* w[(N - k + 1):N]
    lmul!(β, view(w, (N - k + 1):N))
    @views w[1:(N - k)] .= α .* L.coeffs .* v[(k + 1):N] .+ β .* w[1:(N - k)]
    return w
end

# ─── mul! without precomputed coefficients ───────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyPowerOperator{T, false}, v::AbstractVecOrMat) where {T}
    N, k = L.N, L.k

    # @inbounds for i in 1:(N - k)
    #     w[i] = _power_coeff(T, i, k) * v[i + k]
    # end
    # @inbounds for i in (N - k + 1):N
    #     w[i] = zero(eltype(w))
    # end

    fill!(w, zero(eltype(w)))
    @views w[1:(N - k)] .= _power_coeff.(T, 1:(N - k), Ref(k)) .* v[(k + 1):N]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyPowerOperator{T, false}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.N, L.k

    # @inbounds for i in 1:(N - k)
    #     w[i] = α * _power_coeff(T, i, k) * v[i + k] + β * w[i]
    # end
    # @inbounds for i in (N - k + 1):N
    #     w[i] = β * w[i]
    # end

    lmul!(β, view(w, (N - k + 1):N))
    @views w[1:(N - k)] .= α .* _power_coeff.(T, 1:(N - k), Ref(k)) .* v[(k + 1):N] .+ β .* w[1:(N - k)]

    return w
end

# ─── Adjoint: (â†)^k ─────────────────────────────────────────────────────────
# (â†)^k |n⟩ = √((n+k)!/n!) |n+k⟩  ⟹  w[1:k] = 0;  w[i+k] = coeff(i,k) · v[i]

const AdjointDestroyPowerOperator{T} = AdjointOperator{T, <:DestroyPowerOperator{T}} where {T}

# ─── Adjoint mul! with precomputed coefficients ─────────────────────────────

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T, true}}, v::AbstractVecOrMat,
    ) where {T}
    N, k = L.L.N, L.L.k
    # @views w[1:k] .= zero(eltype(w))
    fill!(w, zero(eltype(w)))
    @views w[(k + 1):N] .= L.L.coeffs .* v[1:(N - k)]
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T, true}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.L.N, L.L.k

    # @views w[1:k] .= β .* w[1:k] # This doesn't work on Reactant
    lmul!(β, view(w, 1:k))
    @views w[(k + 1):N] .= α .* L.L.coeffs .* v[1:(N - k)] .+ β .* w[(k + 1):N]
    return w
end

# ─── Adjoint mul! without precomputed coefficients ──────────────────────────

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T, false}}, v::AbstractVecOrMat,
    ) where {T}
    N, k = L.L.N, L.L.k

    # @inbounds for i in 1:k
    #     w[i] = zero(eltype(w))
    # end
    # @inbounds for i in 1:(N - k)
    #     w[i + k] = _power_coeff(T, i, k) * v[i]
    # end

    fill!(w, zero(eltype(w)))
    @views w[(k + 1):N] .= _power_coeff.(T, 1:(N - k), Ref(k)) .* v[1:(N - k)]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T, false}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.L.N, L.L.k

    # @inbounds for i in 1:k
    #     w[i] = β * w[i]
    # end
    # @inbounds for i in 1:(N - k)
    #     w[i + k] = α * _power_coeff(T, i, k) * v[i] + β * w[i + k]
    # end

    lmul!(β, view(w, 1:k))
    @views w[(k + 1):N] .= α .* _power_coeff.(T, 1:(N - k), Ref(k)) .* v[1:(N - k)] .+ β .* w[(k + 1):N]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NormalOrderedOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    NormalOrderedOperator{T,Precomp} <: BosonicOperator{T}

Matrix-free normal-ordered bosonic operator ``(â†)^k â^n``.

Action on Fock basis: ``(â†)^k â^n |m⟩ = \\sqrt{\\frac{m!}{(m-n)!} \\cdot \\frac{(m-n+k)!}{(m-n)!}} |m-n+k⟩``
for ``m ≥ n``, else ``0``.

This is a band operator on diagonal ``k - n``.

When `Precomp = true`, the band coefficients are precomputed and stored.
When `Precomp = false`, they are computed on the fly.
"""
struct NormalOrderedOperator{T, Precomp, VT} <: BosonicOperator{T}
    N::Int
    k::Int   # creation power (â†)^k
    n::Int   # annihilation power â^n
    coeffs::VT

    function NormalOrderedOperator{T, true}(N::Int, k::Int, n::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Creation power k must be positive, got $k"))
        n > 0 || throw(ArgumentError("Annihilation power n must be positive, got $n"))
        len = N - max(k, n)
        len > 0 || throw(ArgumentError("N=$N is too small for (â†)^$k â^$n"))
        coeffs = [_power_coeff(T, i, n) * _power_coeff(T, i, k) for i in 1:len]
        return new{T, true, typeof(coeffs)}(N, k, n, coeffs)
    end

    function NormalOrderedOperator{T, false}(N::Int, k::Int, n::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Creation power k must be positive, got $k"))
        n > 0 || throw(ArgumentError("Annihilation power n must be positive, got $n"))
        len = N - max(k, n)
        len > 0 || throw(ArgumentError("N=$N is too small for (â†)^$k â^$n"))
        return new{T, false, Nothing}(N, k, n, nothing)
    end

    function NormalOrderedOperator(N::Int, k::Int, n::Int, coeffs::Union{Nothing, AbstractVector})
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Creation power k must be positive, got $k"))
        n > 0 || throw(ArgumentError("Annihilation power n must be positive, got $n"))
        if isnothing(coeffs)
            return new{Float64, false, Nothing}(N, k, n, nothing)
        else
            return new{eltype(coeffs), true, typeof(coeffs)}(N, k, n, coeffs)
        end
    end
end

NormalOrderedOperator{T}(N::Int, k::Int, n::Int) where {T} = NormalOrderedOperator{T, true}(N, k, n)
NormalOrderedOperator(N::Int, k::Int, n::Int) = NormalOrderedOperator{Float64}(N, k, n)

Base.size(L::NormalOrderedOperator) = (L.N, L.N)
Base.size(L::NormalOrderedOperator, _::Int) = L.N

islinear(::NormalOrderedOperator) = true
has_adjoint(::NormalOrderedOperator) = true

# adjoint((â†)^k â^n) = (â†)^n â^k  →  swap k and n
function Base.adjoint(L::NormalOrderedOperator{T, Precomp}) where {T, Precomp}
    return NormalOrderedOperator{T, Precomp}(L.N, L.n, L.k)
end

# ─── Helper: coefficient for normal-ordered product ──────────────────────────
# c_j = _power_coeff(T, j, n) * _power_coeff(T, j, k)
# = sqrt(j*(j+1)*...*(j+n-1)) * sqrt(j*(j+1)*...*(j+k-1))

function _normal_ordered_coeff(::Type{T}, j::Int, n::Int, k::Int) where {T}
    return _power_coeff(T, j, n) * _power_coeff(T, j, k)
end

# ─── mul! with precomputed coefficients ──────────────────────────────────────
# (â†)^k â^n maps v[j+n] → w[j+k] with coefficient coeffs[j], for j = 1..N-max(k,n)

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NormalOrderedOperator{T, true}, v::AbstractVecOrMat) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    fill!(w, zero(eltype(w)))
    @views w[1+k:len+k] .= L.coeffs .* v[1+n:len+n]
    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NormalOrderedOperator{T, true}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)

    lmul!(β, w)
    @views w[1+k:len+k] .+= α .* L.coeffs .* v[1+n:len+n]
    return w
end

# ─── mul! without precomputed coefficients ───────────────────────────────────

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NormalOrderedOperator{T, false}, v::AbstractVecOrMat) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    fill!(w, zero(eltype(w)))

    # @inbounds @simd for j in 1:len
    #     w[j + k] = _normal_ordered_coeff(T, j, n, k) * v[j + n]
    # end

    @views w[1+k:len+k] .= _normal_ordered_coeff.(T, 1:len, Ref(n), Ref(k)) .* v[1+n:len+n]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NormalOrderedOperator{T, false}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)

    # @inbounds @simd for i in 1:N
    #     w[i] = β * w[i]
    # end
    # @inbounds @simd for j in 1:len
    #     w[j + k] = α * _normal_ordered_coeff(T, j, n, k) * v[j + n] + w[j + k]
    # end

    lmul!(β, w)
    @views w[1+k:len+k] .+= α .* _normal_ordered_coeff.(T, 1:len, Ref(n), Ref(k)) .* v[1+n:len+n]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Algebraic Simplifications
# ══════════════════════════════════════════════════════════════════════════════

const BosonicOrAdjoint{T} = Union{BosonicOperator{T}, AdjointOperator{T, <:BosonicOperator{T}}} where {T}

# ─── ScaledOperator unwrap rules ─────────────────────────────────────────────
# Peel the scalar so that the inner operator hits the algebraic simplification
# rules below. E.g. `(Δ * a') * a` → `Δ * (a' * a)` → `Δ * NumberOperator`.

function Base.:*(A::ScaledOperator{<:Any, <:Any, <:BosonicOrAdjoint}, B::BosonicOrAdjoint)
    return A.λ * (A.L * B)
end

function Base.:*(A::BosonicOrAdjoint, B::ScaledOperator{<:Any, <:Any, <:BosonicOrAdjoint})
    return (A * B.L) * B.λ
end

# ─── a' * a  →  NumberOperator ───────────────────────────────────────────────
function Base.:*(A::AdjointDestroyOperator{TA, Precomp}, B::DestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.N "Dimension mismatch: a'($(A.L.N)) * a($(B.N))"
    return NumberOperator{promote_type(TA, TB), Precomp}(B.N)
end

# ─── a * a'  →  NumberOperator(shift=1) ──────────────────────────────────────
function Base.:*(A::DestroyOperator{TA, Precomp}, B::AdjointDestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.N == B.L.N "Dimension mismatch: a($(A.N)) * a'($(B.L.N))"
    return NumberOperator{promote_type(TA, TB), Precomp}(A.N; shift = 1)
end

# ─── Destroy × Destroy compositions → DestroyPowerOperator ──────────────────
function Base.:*(A::DestroyOperator{TA, Precomp}, B::DestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.N, 2)
end

function Base.:*(A::DestroyPowerOperator{TA, Precomp}, B::DestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.N, A.k + 1)
end

function Base.:*(A::DestroyOperator{TA, Precomp}, B::DestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.N, 1 + B.k)
end

function Base.:*(A::DestroyPowerOperator{TA, Precomp}, B::DestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.N, A.k + B.k)
end

# ─── Create × Create compositions → adjoint(DestroyPowerOperator) ────────────
function Base.:*(A::AdjointDestroyOperator{TA, Precomp}, B::AdjointDestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.L.N, 2))
end

function Base.:*(A::AdjointDestroyPowerOperator{TA, Precomp}, B::AdjointDestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.L.N, A.L.k + 1))
end

function Base.:*(A::AdjointDestroyOperator{TA, Precomp}, B::AdjointDestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.L.N, 1 + B.L.k))
end

function Base.:*(A::AdjointDestroyPowerOperator{TA, Precomp}, B::AdjointDestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB), Precomp}(A.L.N, A.L.k + B.L.k))
end

# ─── Create × Destroy mixed compositions → NormalOrderedOperator ─────────────
function Base.:*(A::AdjointOperator{TA, <:DestroyPowerOperator{TA, Precomp}}, B::DestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.N "Dimension mismatch: (a')^$(A.L.k)($(A.L.N)) * a^$(B.k)($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB), Precomp}(B.N, A.L.k, B.k)
end

function Base.:*(A::AdjointOperator{TA, <:DestroyPowerOperator{TA, Precomp}}, B::DestroyOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.N "Dimension mismatch: (a')^$(A.L.k)($(A.L.N)) * a($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB), Precomp}(B.N, A.L.k, 1)
end

function Base.:*(A::AdjointOperator{TA, <:DestroyOperator{TA, Precomp}}, B::DestroyPowerOperator{TB, Precomp}) where {TA, TB, Precomp}
    @assert A.L.N == B.N "Dimension mismatch: a'($(A.L.N)) * a^$(B.k)($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB), Precomp}(B.N, 1, B.k)
end

# ─── Powers ──────────────────────────────────────────────────────────────────
function Base.:^(a::DestroyOperator{T, Precomp}, k::Integer) where {T, Precomp}
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return DestroyPowerOperator{T, Precomp}(a.N, k)
end

function Base.:^(a::AdjointDestroyOperator{T, Precomp}, k::Integer) where {T, Precomp}
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return adjoint(DestroyPowerOperator{T, Precomp}(a.L.N, k))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  concretize: convert to sparse matrix
# ═══════════════════════════════════════════════════════════════════════════════

function concretize(L::DestroyOperator{T, true}) where {T}
    return spdiagm(1 => L.coeffs)
end

function concretize(L::DestroyOperator{T, false}) where {T}
    N = L.N
    return spdiagm(1 => [_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::AdjointOperator{T, <:DestroyOperator{T, true}}) where {T}
    return spdiagm(-1 => L.L.coeffs)
end

function concretize(L::AdjointOperator{T, <:DestroyOperator{T, false}}) where {T}
    N = L.L.N
    return spdiagm(-1 => [_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::NumberOperator{T, true}) where {T}
    return spdiagm(0 => L.coeffs)
end

function concretize(L::NumberOperator{T, false}) where {T}
    N = L.N
    return spdiagm(0 => [real(T)(i - 1 + L.shift) for i in 1:N])
end

function concretize(L::DestroyPowerOperator{T, true}) where {T}
    return spdiagm(L.k => L.coeffs)
end

function concretize(L::DestroyPowerOperator{T, false}) where {T}
    N, k = L.N, L.k
    return spdiagm(k => [_power_coeff(T, i, k) for i in 1:(N - k)])
end

function concretize(L::AdjointOperator{T, <:DestroyPowerOperator{T, true}}) where {T}
    return spdiagm(-L.L.k => L.L.coeffs)
end

function concretize(L::AdjointOperator{T, <:DestroyPowerOperator{T, false}}) where {T}
    N, k = L.L.N, L.L.k
    return spdiagm(-k => [_power_coeff(T, i, k) for i in 1:(N - k)])
end

function concretize(L::NormalOrderedOperator{T, true}) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    rows = (1:len) .+ k
    cols = (1:len) .+ n
    return sparse(rows, cols, L.coeffs, N, N)
end

function concretize(L::NormalOrderedOperator{T, false}) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    rows = (1:len) .+ k
    cols = (1:len) .+ n
    coeffs = [_normal_ordered_coeff(T, j, n, k) for j in 1:len]
    return sparse(rows, cols, coeffs, N, N)
end
