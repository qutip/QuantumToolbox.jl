# ──────────────────────────────────────────────────────────────────────────────
# Bosonic matrix-free quantum operators
# ──────────────────────────────────────────────────────────────────────────────

# ─── Type-aware coefficient helpers ──────────────────────────────────────────
# Avoid sqrt(::Int) → Float64 promotion when working with Float32 / other types.

_sqrt_coeff(::Type{T}, i::Int) where {T} = sqrt(T(i))
_sqrt_coeff(::Type{Complex{T}}, i::Int) where {T} = _sqrt_coeff(T, i)

function _power_coeff(::Type{T}, i::Int, k::Int) where {T}
    c = one(T)
    for m in i:(i + k - 1)
        c *= T(m)
    end
    return sqrt(c)
end
_power_coeff(::Type{Complex{T}}, i::Int, k::Int) where {T} = _power_coeff(T, i, k)

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
struct DestroyOperator{T} <: AbstractSciMLOperator{T}
    N::Int

    function DestroyOperator{T}(N::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T}(N)
    end
end

DestroyOperator(N::Int) = DestroyOperator{Float64}(N)

Base.size(L::DestroyOperator) = (L.N, L.N)
Base.size(L::DestroyOperator, n::Int) = L.N

islinear(::DestroyOperator) = true
has_adjoint(::DestroyOperator) = true

# â |n⟩ = √n |n-1⟩  ⟹  w[i] = √i · v[i+1]  for i = 1…N-1;  w[N] = 0
function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyOperator{T}, v::AbstractVecOrMat) where {T}
    N = L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - 1)
        w[i] = _sqrt_coeff(T, i) * v[i + 1]
    end
    @inbounds w[N] = zero(eltype(w))
    return w
end

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::DestroyOperator{T}, v::AbstractVecOrMat, α, β,
) where {T}
    N = L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - 1)
        w[i] = α * _sqrt_coeff(T, i) * v[i + 1] + β * w[i]
    end
    @inbounds w[N] = β * w[N]
    return w
end

# ─── Adjoint: creation operator â† ───────────────────────────────────────────
# â† |n⟩ = √(n+1) |n+1⟩  ⟹  w[1] = 0;  w[i+1] = √i · v[i]  for i = 1…N-1

const AdjointDestroyOperator{T} = AdjointOperator{T, DestroyOperator{T}} where T

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyOperator{T}, v::AbstractVecOrMat,
) where {T}
    N = L.L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds w[1] = zero(eltype(w))
    @inbounds for i in 1:(N - 1)
        w[i + 1] = _sqrt_coeff(T, i) * v[i]
    end
    return w
end

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyOperator{T}, v::AbstractVecOrMat, α, β,
) where {T}
    N = L.L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds w[1] = β * w[1]
    @inbounds for i in 1:(N - 1)
        w[i + 1] = α * _sqrt_coeff(T, i) * v[i] + β * w[i + 1]
    end
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
struct NumberOperator{T} <: AbstractSciMLOperator{T}
    N::Int
    shift::Int

    function NumberOperator{T}(N::Int; shift::Int = 0) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new{T}(N, shift)
    end
end

NumberOperator(N::Int; shift::Int = 0) = NumberOperator{Float64}(N; shift)

Base.size(L::NumberOperator) = (L.N, L.N)
Base.size(L::NumberOperator, n::Int) = L.N

islinear(::NumberOperator) = true
has_adjoint(::NumberOperator) = true

# NumberOperator is Hermitian: adjoint returns itself
Base.adjoint(L::NumberOperator) = L

function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NumberOperator{T}, v::AbstractVecOrMat) where {T}
    N, shift = L.N, L.shift
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:N
        w[i] = T(i - 1 + shift) * v[i]
    end
    return w
end

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::NumberOperator{T}, v::AbstractVecOrMat, α, β,
) where {T}
    N, shift = L.N, L.shift
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:N
        w[i] = α * T(i - 1 + shift) * v[i] + β * w[i]
    end
    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyPowerOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyPowerOperator{T} <: AbstractSciMLOperator{T}

Matrix-free ``â^k`` (k-th power of the bosonic annihilation operator).

Action on Fock basis: ``â^k |n⟩ = \\sqrt{n! / (n-k)!}\\, |n-k⟩`` for ``n ≥ k``, else ``0``.

Computed in a single O(N) pass.
"""
struct DestroyPowerOperator{T} <: AbstractSciMLOperator{T}
    N::Int
    k::Int

    function DestroyPowerOperator{T}(N::Int, k::Int) where {T}
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        return new{T}(N, k)
    end
end

DestroyPowerOperator(N::Int, k::Int) = DestroyPowerOperator{Float64}(N, k)

Base.size(L::DestroyPowerOperator) = (L.N, L.N)
Base.size(L::DestroyPowerOperator, n::Int) = L.N

islinear(::DestroyPowerOperator) = true
has_adjoint(::DestroyPowerOperator) = true

# â^k |n⟩ = √(n!/(n-k)!) |n-k⟩  ⟹  w[i] = coeff(i,k) · v[i+k]
function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyPowerOperator{T}, v::AbstractVecOrMat) where {T}
    N, k = L.N, L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - k)
        w[i] = _power_coeff(T, i, k) * v[i + k]
    end
    @inbounds for i in (N - k + 1):N
        w[i] = zero(eltype(w))
    end
    return w
end

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::DestroyPowerOperator{T}, v::AbstractVecOrMat, α, β,
) where {T}
    N, k = L.N, L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - k)
        w[i] = α * _power_coeff(T, i, k) * v[i + k] + β * w[i]
    end
    @inbounds for i in (N - k + 1):N
        w[i] = β * w[i]
    end
    return w
end

# ─── Adjoint: (â†)^k ─────────────────────────────────────────────────────────
# (â†)^k |n⟩ = √((n+k)!/n!) |n+k⟩  ⟹  w[1:k] = 0;  w[i+k] = coeff(i,k) · v[i]

const AdjointDestroyPowerOperator{T} = AdjointOperator{T, DestroyPowerOperator{T}} where T

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyPowerOperator{T}, v::AbstractVecOrMat,
) where {T}
    N, k = L.L.N, L.L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:k
        w[i] = zero(eltype(w))
    end
    @inbounds for i in 1:(N - k)
        w[i + k] = _power_coeff(T, i, k) * v[i]
    end
    return w
end

function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyPowerOperator{T}, v::AbstractVecOrMat, α, β,
) where {T}
    N, k = L.L.N, L.L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:k
        w[i] = β * w[i]
    end
    @inbounds for i in 1:(N - k)
        w[i + k] = α * _power_coeff(T, i, k) * v[i] + β * w[i + k]
    end
    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Algebraic Simplifications
# ══════════════════════════════════════════════════════════════════════════════

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
    return spdiagm(1 => T[_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::AdjointDestroyOperator{T}) where {T}
    N = L.L.N
    return spdiagm(-1 => T[_sqrt_coeff(T, i) for i in 1:(N - 1)])
end

function concretize(L::NumberOperator{T}) where {T}
    N = L.N
    return spdiagm(0 => T[T(i - 1 + L.shift) for i in 1:N])
end

function concretize(L::DestroyPowerOperator{T}) where {T}
    N, k = L.N, L.k
    return spdiagm(k => T[_power_coeff(T, i, k) for i in 1:(N - k)])
end

function concretize(L::AdjointDestroyPowerOperator{T}) where {T}
    N, k = L.L.N, L.L.k
    return spdiagm(-k => T[_power_coeff(T, i, k) for i in 1:(N - k)])
end
