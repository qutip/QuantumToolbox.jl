# ──────────────────────────────────────────────────────────────────────────────
# Bosonic matrix-free quantum operators
# ──────────────────────────────────────────────────────────────────────────────

"""
    BosonicOperator{T} <: AbstractSciMLOperator{T}

Abstract supertype for all matrix-free bosonic operators.
"""
abstract type BosonicOperator{T} <: AbstractSciMLOperator{T} end

const BosonicOrAdjoint{T} = Union{BosonicOperator{T}, AdjointOperator{T, <:BosonicOperator{T}}} where {T}

# ─── Type-aware coefficient helpers ──────────────────────────────────────────
# Avoid sqrt(::Int) → Float64 promotion when working with Float32 / other types.

_sqrt_coeff(::Type{T}, i::Integer) where {T} = sqrt(real(T)(i))

function _power_coeff(::Type{T}, i::Integer, k::Integer) where {T}
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
    DestroyOperator{T} <: AbstractSciMLOperator{T}

Matrix-free bosonic annihilation (destroy) operator.

Action on a Fock-basis vector: ``â |n⟩ = √n |n-1⟩``.

The creation operator ``â†`` is obtained via `adjoint(a)`, which returns
an `AdjointOperator{T, DestroyOperator}` — no separate type is needed.
"""
struct DestroyOperator{T} <: BosonicOperator{T}
    N::Int

    function DestroyOperator{T}(N::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T}(N)
    end
end

DestroyOperator(N::Int) = DestroyOperator{Float64}(N)

Base.size(L::DestroyOperator) = (L.N, L.N)
Base.size(L::DestroyOperator, n::Int) = size(L)[n]

islinear(::DestroyOperator) = true
has_adjoint(::DestroyOperator) = true

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyOperator{T}, v::AbstractVecOrMat) where {T}
    N = L.N

    fill!(w, zero(eltype(w)))
    @views w[1:(N - 1)] .= _sqrt_coeff.(T, 1:(N - 1)) .* v[2:N]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyOperator{T}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.N

    lmul!(β, w)
    @views w[1:(N - 1)] .= α .* _sqrt_coeff.(T, 1:(N - 1)) .* v[2:N] .+ β .* w[1:(N - 1)]

    return w
end

# ─── Adjoint: creation operator â† ───────────────────────────────────────────
# â† |n⟩ = √(n+1) |n+1⟩  ⟹  w[1] = 0;  w[i+1] = √i · v[i]  for i = 1…N-1

const AdjointDestroyOperator{T} = AdjointOperator{T, <:DestroyOperator{T}} where {T}

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T}}, v::AbstractVecOrMat,
    ) where {T}
    N = L.L.N

    fill!(w, zero(eltype(w)))
    @views w[2:N] .= _sqrt_coeff.(T, 1:(N - 1)) .* v[1:(N - 1)]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyOperator{T}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N = L.L.N

    lmul!(β, w)
    @views w[2:N] .= α .* _sqrt_coeff.(T, 1:(N - 1)) .* v[1:(N - 1)] .+ β .* w[2:N]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NumberOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    NumberOperator{T} <: AbstractSciMLOperator{T}

Matrix-free bosonic number operator with optional shift.

Action: ``w_i = (i - 1 + \\text{shift}) \\cdot v_i``.

    - `shift = 0` represents ``\\hat{n} = â†â`` with eigenvalues ``0, 1, 2, …``
  - `shift = 1` represents ``ââ†`` with eigenvalues ``1, 2, 3, …``
"""
struct NumberOperator{T} <: BosonicOperator{T}
    N::Int
    shift::T

    function NumberOperator{T1}(N::Int; shift::T2 = zero(T1)) where {T1, T2}
        T2 <: T1 || throw(ArgumentError("Shift type $T2 is not a subtype of $T1"))
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T1}(N, shift)
    end
end

NumberOperator(N::Int; shift::T = zero(Float64)) where {T} = NumberOperator{T}(N; shift)

Base.size(L::NumberOperator) = (L.N, L.N)
Base.size(L::NumberOperator, n::Int) = size(L)[n]

islinear(::NumberOperator) = true
has_adjoint(::NumberOperator) = true

# NumberOperator is Hermitian: adjoint returns itself
Base.adjoint(L::NumberOperator) = L

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NumberOperator{T}, v::AbstractVecOrMat) where {T}
    N, shift = L.N, L.shift

    fill!(w, zero(eltype(w)))
    w .= (real(T).(0:(N - 1)) .+ shift) .* v

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NumberOperator{T}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, shift = L.N, L.shift

    lmul!(β, w)
    w .= α .* (real(T).(0:(N - 1)) .+ shift) .* v .+ β .* w

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyPowerOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyPowerOperator{T} <: AbstractSciMLOperator{T}

Matrix-free ``â^k`` (k-th power of the bosonic annihilation operator).

Action on Fock basis: ``â^k |n⟩ = \\sqrt{n! / (n-k)!}\\, |n-k⟩`` for ``n ≥ k``, else ``0``.
"""
struct DestroyPowerOperator{T} <: BosonicOperator{T}
    N::Int
    k::Int32

    function DestroyPowerOperator{T}(N::Int, k::Integer) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        return new{T}(N, Int32(k))
    end
end

DestroyPowerOperator(N::Int, k::Integer) = DestroyPowerOperator{Float64}(N, k)

Base.size(L::DestroyPowerOperator) = (L.N, L.N)
Base.size(L::DestroyPowerOperator, n::Int) = size(L)[n]

islinear(::DestroyPowerOperator) = true
has_adjoint(::DestroyPowerOperator) = true

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyPowerOperator{T}, v::AbstractVecOrMat) where {T}
    N, k = L.N, L.k

    fill!(w, zero(eltype(w)))
    @views w[1:(N - k)] .= _power_coeff.(T, 1:(N - k), Ref(k)) .* v[(k + 1):N]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::DestroyPowerOperator{T}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.N, L.k

    lmul!(β, view(w, (N - k + 1):N))
    @views w[1:(N - k)] .= α .* _power_coeff.(T, 1:(N - k), Ref(k)) .* v[(k + 1):N] .+ β .* w[1:(N - k)]

    return w
end

# ─── Adjoint: (â†)^k ─────────────────────────────────────────────────────────
# (â†)^k |n⟩ = √((n+k)!/n!) |n+k⟩  ⟹  w[1:k] = 0;  w[i+k] = coeff(i,k) · v[i]

const AdjointDestroyPowerOperator{T} = AdjointOperator{T, <:DestroyPowerOperator{T}} where {T}

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T}}, v::AbstractVecOrMat,
    ) where {T}
    N, k = L.L.N, L.L.k

    fill!(w, zero(eltype(w)))
    @views w[(k + 1):N] .= _power_coeff.(T, 1:(N - k), Ref(k)) .* v[1:(N - k)]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::AdjointOperator{T, <:DestroyPowerOperator{T}}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k = L.L.N, L.L.k

    lmul!(β, view(w, 1:k))
    @views w[(k + 1):N] .= α .* _power_coeff.(T, 1:(N - k), Ref(k)) .* v[1:(N - k)] .+ β .* w[(k + 1):N]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NormalOrderedOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    NormalOrderedOperator{T} <: BosonicOperator{T}

Matrix-free normal-ordered bosonic operator ``(â†)^k â^n``.

Action on Fock basis: ``(â†)^k â^n |m⟩ = \\sqrt{\\frac{m!}{(m-n)!} \\cdot \\frac{(m-n+k)!}{(m-n)!}} |m-n+k⟩``
for ``m ≥ n``, else ``0``.

This is a band operator on diagonal ``k - n``.
"""
struct NormalOrderedOperator{T} <: BosonicOperator{T}
    N::Int
    k::Int32
    n::Int32

    function NormalOrderedOperator{T}(N::Int, k::Integer, n::Integer) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Creation power k must be positive, got $k"))
        n > 0 || throw(ArgumentError("Annihilation power n must be positive, got $n"))
        len = N - max(k, n)
        len > 0 || throw(ArgumentError("N=$N is too small for (â†)^$k â^$n"))
        return new{T}(N, Int32(k), Int32(n))
    end
end

NormalOrderedOperator(N::Int, k::Integer, n::Integer) = NormalOrderedOperator{Float64}(N, k, n)

Base.size(L::NormalOrderedOperator) = (L.N, L.N)
Base.size(L::NormalOrderedOperator, _::Int) = size(L)[1]

islinear(::NormalOrderedOperator) = true
has_adjoint(::NormalOrderedOperator) = true

# adjoint((â†)^k â^n) = (â†)^n â^k  →  swap k and n
function Base.adjoint(L::NormalOrderedOperator{T}) where {T}
    return NormalOrderedOperator{T}(L.N, L.n, L.k)
end

# ─── Helper: coefficient for normal-ordered product ──────────────────────────
# c_j = _power_coeff(T, j, n) * _power_coeff(T, j, k)
# = sqrt(j*(j+1)*...*(j+n-1)) * sqrt(j*(j+1)*...*(j+k-1))

function _normal_ordered_coeff(::Type{T}, j::Int, n::Integer, k::Integer) where {T}
    return _power_coeff(T, j, n) * _power_coeff(T, j, k)
end

# (â†)^k â^n maps v[j+n] → w[j+k] with coefficient coeff(j), for j = 1..N-max(k,n)

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NormalOrderedOperator{T}, v::AbstractVecOrMat) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    fill!(w, zero(eltype(w)))

    @views w[(k + 1):(len + k)] .= _normal_ordered_coeff.(T, 1:len, Ref(n), Ref(k)) .* v[(n + 1):(len + n)]

    return w
end

function LinearAlgebra.mul!(
        w::AbstractVecOrMat, L::NormalOrderedOperator{T}, v::AbstractVecOrMat, α, β,
    ) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)

    lmul!(β, w)
    @views w[(k + 1):(len + k)] .+= α .* _normal_ordered_coeff.(T, 1:len, Ref(n), Ref(k)) .* v[(n + 1):(len + n)]

    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Algebraic Simplifications
# ══════════════════════════════════════════════════════════════════════════════

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
function Base.:*(A::AdjointDestroyOperator{TA}, B::DestroyOperator{TB}) where {TA, TB}
    @assert A.L.N == B.N "Dimension mismatch: a'($(A.L.N)) * a($(B.N))"
    return NumberOperator{promote_type(TA, TB)}(B.N)
end

# ─── a * a'  →  NumberOperator(shift=1) ──────────────────────────────────────
function Base.:*(A::DestroyOperator{TA}, B::AdjointDestroyOperator{TB}) where {TA, TB}
    @assert A.N == B.L.N "Dimension mismatch: a($(A.N)) * a'($(B.L.N))"
    return NumberOperator{promote_type(TA, TB)}(A.N; shift = 1)
end

# ─── Destroy × Destroy compositions → DestroyPowerOperator ──────────────────
function Base.:*(A::DestroyOperator{TA}, B::DestroyOperator{TB}) where {TA, TB}
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB)}(A.N, 2)
end

function Base.:*(A::DestroyPowerOperator{TA}, B::DestroyOperator{TB}) where {TA, TB}
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB)}(A.N, A.k + 1)
end

function Base.:*(A::DestroyOperator{TA}, B::DestroyPowerOperator{TB}) where {TA, TB}
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB)}(A.N, 1 + B.k)
end

function Base.:*(A::DestroyPowerOperator{TA}, B::DestroyPowerOperator{TB}) where {TA, TB}
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator{promote_type(TA, TB)}(A.N, A.k + B.k)
end

# ─── Create × Create compositions → adjoint(DestroyPowerOperator) ────────────
function Base.:*(A::AdjointDestroyOperator{TA}, B::AdjointDestroyOperator{TB}) where {TA, TB}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB)}(A.L.N, 2))
end

function Base.:*(A::AdjointDestroyPowerOperator{TA}, B::AdjointDestroyOperator{TB}) where {TA, TB}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB)}(A.L.N, A.L.k + 1))
end

function Base.:*(A::AdjointDestroyOperator{TA}, B::AdjointDestroyPowerOperator{TB}) where {TA, TB}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB)}(A.L.N, 1 + B.L.k))
end

function Base.:*(A::AdjointDestroyPowerOperator{TA}, B::AdjointDestroyPowerOperator{TB}) where {TA, TB}
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator{promote_type(TA, TB)}(A.L.N, A.L.k + B.L.k))
end

# ─── Create × Destroy mixed compositions → NormalOrderedOperator ─────────────
function Base.:*(A::AdjointOperator{TA, <:DestroyPowerOperator{TA}}, B::DestroyPowerOperator{TB}) where {TA, TB}
    @assert A.L.N == B.N "Dimension mismatch: (a')^$(A.L.k)($(A.L.N)) * a^$(B.k)($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB)}(B.N, A.L.k, B.k)
end

function Base.:*(A::AdjointOperator{TA, <:DestroyPowerOperator{TA}}, B::DestroyOperator{TB}) where {TA, TB}
    @assert A.L.N == B.N "Dimension mismatch: (a')^$(A.L.k)($(A.L.N)) * a($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB)}(B.N, A.L.k, 1)
end

function Base.:*(A::AdjointOperator{TA, <:DestroyOperator{TA}}, B::DestroyPowerOperator{TB}) where {TA, TB}
    @assert A.L.N == B.N "Dimension mismatch: a'($(A.L.N)) * a^$(B.k)($(B.N))"
    return NormalOrderedOperator{promote_type(TA, TB)}(B.N, 1, B.k)
end

# ─── Powers ──────────────────────────────────────────────────────────────────
function Base.:^(a::DestroyOperator{T}, k::Integer) where {T}
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return DestroyPowerOperator{T}(a.N, k)
end

function Base.:^(a::AdjointDestroyOperator{T}, k::Integer) where {T}
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return adjoint(DestroyPowerOperator{T}(a.L.N, k))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  concretize: convert to sparse matrix
# ═══════════════════════════════════════════════════════════════════════════════

function concretize(L::DestroyOperator{T}) where {T}
    N = L.N
    return spdiagm(1 => [_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::AdjointOperator{T, <:DestroyOperator{T}}) where {T}
    N = L.L.N
    return spdiagm(-1 => [_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::NumberOperator{T}) where {T}
    N = L.N
    return spdiagm(0 => [real(T)(i - 1 + L.shift) for i in 1:N])
end

function concretize(L::DestroyPowerOperator{T}) where {T}
    N, k = L.N, L.k
    return spdiagm(k => [_power_coeff(T, i, k) for i in 1:(N - k)])
end

function concretize(L::AdjointOperator{T, <:DestroyPowerOperator{T}}) where {T}
    N, k = L.L.N, L.L.k
    return spdiagm(-k => [_power_coeff(T, i, k) for i in 1:(N - k)])
end

function concretize(L::NormalOrderedOperator{T}) where {T}
    N, k, n = L.N, L.k, L.n
    len = N - max(k, n)
    rows = (1:len) .+ k
    cols = (1:len) .+ n
    coeffs = [_normal_ordered_coeff(T, j, n, k) for j in 1:len]
    return sparse(rows, cols, coeffs, N, N)
end
