# ──────────────────────────────────────────────────────────────────────────────
# Bosonic matrix-free quantum operators
# ──────────────────────────────────────────────────────────────────────────────

# ─── Type-aware coefficient helpers ──────────────────────────────────────────
# Avoid sqrt(::Int) → Float64 promotion when working with Float32 / other types.

@inline _sqrt_coeff(::Type{T}, i::Int) where {T <: Real} = sqrt(T(i))
@inline _sqrt_coeff(v::AbstractVecOrMat, i::Int) = _sqrt_coeff(real(eltype(v)), i)

@inline function _power_coeff(::Type{T}, i::Int, k::Int) where {T <: Real}
    c = one(T)
    for m in i:(i + k - 1)
        c *= T(m)
    end
    return sqrt(c)
end
@inline _power_coeff(v::AbstractVecOrMat, i::Int, k::Int) = _power_coeff(real(eltype(v)), i, k)

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyOperator <: AbstractSciMLOperator{Bool}

Matrix-free bosonic annihilation (destroy) operator.

Action on a Fock-basis vector: ``â |n⟩ = √n |n-1⟩``.

The creation operator ``â†`` is obtained via `adjoint(a)`, which returns
an `AdjointOperator{Bool, DestroyOperator}` — no separate type is needed.
"""
struct DestroyOperator <: AbstractSciMLOperator{Bool}
    N::Int

    function DestroyOperator(N::Int)
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new(N)
    end
end

Base.size(L::DestroyOperator) = (L.N, L.N)
Base.size(L::DestroyOperator, n::Int) = L.N

islinear(::DestroyOperator) = true
has_adjoint(::DestroyOperator) = true

# â |n⟩ = √n |n-1⟩  ⟹  w[i] = √i · v[i+1]  for i = 1…N-1;  w[N] = 0
@inline function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyOperator, v::AbstractVecOrMat)
    N = L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - 1)
        w[i] = _sqrt_coeff(v, i) * v[i + 1]
    end
    @inbounds w[N] = zero(eltype(w))
    return w
end

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::DestroyOperator, v::AbstractVecOrMat, α, β,
)
    N = L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - 1)
        w[i] = α * _sqrt_coeff(v, i) * v[i + 1] + β * w[i]
    end
    @inbounds w[N] = β * w[N]
    return w
end

# ─── Adjoint: creation operator â† ───────────────────────────────────────────
# â† |n⟩ = √(n+1) |n+1⟩  ⟹  w[1] = 0;  w[i+1] = √i · v[i]  for i = 1…N-1

const AdjointDestroyOperator{T} = AdjointOperator{T, DestroyOperator}

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyOperator, v::AbstractVecOrMat,
)
    N = L.L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds w[1] = zero(eltype(w))
    @inbounds for i in 1:(N - 1)
        w[i + 1] = _sqrt_coeff(v, i) * v[i]
    end
    return w
end

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyOperator, v::AbstractVecOrMat, α, β,
)
    N = L.L.N
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds w[1] = β * w[1]
    @inbounds for i in 1:(N - 1)
        w[i + 1] = α * _sqrt_coeff(v, i) * v[i] + β * w[i + 1]
    end
    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NumberOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    NumberOperator <: AbstractSciMLOperator{Bool}

Matrix-free bosonic number operator with optional shift.

Action: ``w_i = (i - 1 + \\text{shift}) \\cdot v_i``.

  - `shift = 0` represents ``\\hat{n} = â†â`` with eigenvalues ``0, 1, 2, …``
  - `shift = 1` represents ``ââ†`` with eigenvalues ``1, 2, 3, …``
"""
struct NumberOperator <: AbstractSciMLOperator{Bool}
    N::Int
    shift::Int

    function NumberOperator(N::Int; shift::Int = 0)
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        return new(N, shift)
    end
end

Base.size(L::NumberOperator) = (L.N, L.N)
Base.size(L::NumberOperator, n::Int) = L.N

islinear(::NumberOperator) = true
has_adjoint(::NumberOperator) = true

# NumberOperator is Hermitian: adjoint returns itself
Base.adjoint(L::NumberOperator) = L

@inline function LinearAlgebra.mul!(w::AbstractVecOrMat, L::NumberOperator, v::AbstractVecOrMat)
    N, shift = L.N, L.shift
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:N
        w[i] = (i - 1 + shift) * v[i]
    end
    return w
end

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::NumberOperator, v::AbstractVecOrMat, α, β,
)
    N, shift = L.N, L.shift
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:N
        w[i] = α * (i - 1 + shift) * v[i] + β * w[i]
    end
    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  DestroyPowerOperator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    DestroyPowerOperator <: AbstractSciMLOperator{Bool}

Matrix-free ``â^k`` (k-th power of the bosonic annihilation operator).

Action on Fock basis: ``â^k |n⟩ = \\sqrt{n! / (n-k)!}\\, |n-k⟩`` for ``n ≥ k``, else ``0``.

Computed in a single O(N) pass.
"""
struct DestroyPowerOperator <: AbstractSciMLOperator{Bool}
    N::Int
    k::Int

    function DestroyPowerOperator(N::Int, k::Int)
        N > 0 || throw(ArgumentError("Hilbert space dimension N must be positive, got $N"))
        k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
        return new(N, k)
    end
end

Base.size(L::DestroyPowerOperator) = (L.N, L.N)
Base.size(L::DestroyPowerOperator, n::Int) = L.N

islinear(::DestroyPowerOperator) = true
has_adjoint(::DestroyPowerOperator) = true

# â^k |n⟩ = √(n!/(n-k)!) |n-k⟩  ⟹  w[i] = coeff(i,k) · v[i+k]
@inline function LinearAlgebra.mul!(w::AbstractVecOrMat, L::DestroyPowerOperator, v::AbstractVecOrMat)
    N, k = L.N, L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - k)
        w[i] = _power_coeff(v, i, k) * v[i + k]
    end
    @inbounds for i in (N - k + 1):N
        w[i] = zero(eltype(w))
    end
    return w
end

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::DestroyPowerOperator, v::AbstractVecOrMat, α, β,
)
    N, k = L.N, L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:(N - k)
        w[i] = α * _power_coeff(v, i, k) * v[i + k] + β * w[i]
    end
    @inbounds for i in (N - k + 1):N
        w[i] = β * w[i]
    end
    return w
end

# ─── Adjoint: (â†)^k ─────────────────────────────────────────────────────────
# (â†)^k |n⟩ = √((n+k)!/n!) |n+k⟩  ⟹  w[1:k] = 0;  w[i+k] = coeff(i,k) · v[i]

const AdjointDestroyPowerOperator{T} = AdjointOperator{T, DestroyPowerOperator}

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyPowerOperator, v::AbstractVecOrMat,
)
    N, k = L.L.N, L.L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:k
        w[i] = zero(eltype(w))
    end
    @inbounds for i in 1:(N - k)
        w[i + k] = _power_coeff(v, i, k) * v[i]
    end
    return w
end

@inline function LinearAlgebra.mul!(
    w::AbstractVecOrMat, L::AdjointDestroyPowerOperator, v::AbstractVecOrMat, α, β,
)
    N, k = L.L.N, L.L.k
    @assert size(v, 1) == N "Input dimension mismatch: expected $N, got $(size(v, 1))"
    @assert size(w, 1) == N "Output dimension mismatch: expected $N, got $(size(w, 1))"
    @inbounds for i in 1:k
        w[i] = β * w[i]
    end
    @inbounds for i in 1:(N - k)
        w[i + k] = α * _power_coeff(v, i, k) * v[i] + β * w[i + k]
    end
    return w
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Algebraic Simplifications
# ═══════════════════════════════════════════════════════════════════════════════

# ─── a' * a  →  NumberOperator ───────────────────────────────────────────────
function Base.:*(A::AdjointDestroyOperator, B::DestroyOperator)
    @assert A.L.N == B.N "Dimension mismatch: a'($(A.L.N)) * a($(B.N))"
    return NumberOperator(B.N)
end

# ─── a * a'  →  NumberOperator(shift=1) ──────────────────────────────────────
function Base.:*(A::DestroyOperator, B::AdjointDestroyOperator)
    @assert A.N == B.L.N "Dimension mismatch: a($(A.N)) * a'($(B.L.N))"
    return NumberOperator(A.N; shift = 1)
end

# ─── Destroy × Destroy compositions → DestroyPowerOperator ──────────────────
function Base.:*(A::DestroyOperator, B::DestroyOperator)
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a($(B.N))"
    return DestroyPowerOperator(A.N, 2)
end

function Base.:*(A::DestroyPowerOperator, B::DestroyOperator)
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a($(B.N))"
    return DestroyPowerOperator(A.N, A.k + 1)
end

function Base.:*(A::DestroyOperator, B::DestroyPowerOperator)
    @assert A.N == B.N "Dimension mismatch: a($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator(A.N, 1 + B.k)
end

function Base.:*(A::DestroyPowerOperator, B::DestroyPowerOperator)
    @assert A.N == B.N "Dimension mismatch: a^$(A.k)($(A.N)) * a^$(B.k)($(B.N))"
    return DestroyPowerOperator(A.N, A.k + B.k)
end

# ─── Create × Create compositions → adjoint(DestroyPowerOperator) ────────────
function Base.:*(A::AdjointDestroyOperator, B::AdjointDestroyOperator)
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator(A.L.N, 2))
end

function Base.:*(A::AdjointDestroyPowerOperator, B::AdjointDestroyOperator)
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator(A.L.N, A.L.k + 1))
end

function Base.:*(A::AdjointDestroyOperator, B::AdjointDestroyPowerOperator)
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator(A.L.N, 1 + B.L.k))
end

function Base.:*(A::AdjointDestroyPowerOperator, B::AdjointDestroyPowerOperator)
    @assert A.L.N == B.L.N "Dimension mismatch"
    return adjoint(DestroyPowerOperator(A.L.N, A.L.k + B.L.k))
end

# ─── Powers ──────────────────────────────────────────────────────────────────
function Base.:^(a::DestroyOperator, k::Integer)
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return DestroyPowerOperator(a.N, k)
end

function Base.:^(a::AdjointDestroyOperator, k::Integer)
    k > 0 || throw(ArgumentError("Power k must be positive, got $k"))
    k == 1 && return a
    return adjoint(DestroyPowerOperator(a.L.N, k))
end

# ═══════════════════════════════════════════════════════════════════════════════
#  concretize: convert to sparse matrix
# ═══════════════════════════════════════════════════════════════════════════════

function concretize(L::DestroyOperator)
    N = L.N
    return spdiagm(1 => ComplexF64[sqrt(i) for i in 1:(N - 1)])
end

function concretize(L::AdjointDestroyOperator)
    N = L.L.N
    return spdiagm(-1 => ComplexF64[sqrt(i) for i in 1:(N - 1)])
end

function concretize(L::NumberOperator)
    N = L.N
    return spdiagm(0 => ComplexF64[i - 1 + L.shift for i in 1:N])
end

function concretize(L::DestroyPowerOperator)
    N, k = L.N, L.k
    return spdiagm(k => ComplexF64[_power_coeff(Float64, i, k) for i in 1:(N - k)])
end

function concretize(L::AdjointDestroyPowerOperator)
    N, k = L.L.N, L.L.k
    return spdiagm(-k => ComplexF64[_power_coeff(Float64, i, k) for i in 1:(N - k)])
end
