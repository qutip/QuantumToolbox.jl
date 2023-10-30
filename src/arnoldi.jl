export ArnoldiSpace, arnoldi, arnoldi!, arnoldi_init!, arnoldi_step!
export expv!, expv

struct ArnoldiSpace{VT<:AbstractMatrix{<: BlasFloat}, HT<:AbstractMatrix{<: BlasFloat}, mT<:Integer}
    V::VT
    H::HT
    m::mT
end

function Base.copy(AS::ArnoldiSpace{<:AbstractMatrix{T1}, <:AbstractMatrix{T1}}) where T1 <: BlasFloat
    ArnoldiSpace(copy(AS.V), copy(AS.H), AS.m)
end

function Base.deepcopy(AS::ArnoldiSpace{<:AbstractMatrix{T1}, <:AbstractMatrix{T1}}) where T1 <: BlasFloat
    ArnoldiSpace(deepcopy(AS.V), deepcopy(AS.H), AS.m)
end

function arnoldi_init!(A, b::AbstractVector{T}, V::AbstractMatrix{T}, H::AbstractMatrix{T}) where T <: BlasFloat
    v₁ = view(V, :, 1)
    v₂ = view(V, :, 2)
    v₁ .= b
    normalize!(v₁)
    mul!(v₂, A, v₁)
    H[1,1] = dot(v₁, v₂)
    axpy!(-H[1,1], v₁, v₂)
    H[2,1] = norm(v₂)
    @. v₂ = v₂ / H[2,1]
end
  
function arnoldi_step!(A, V::AbstractMatrix{T}, H::AbstractMatrix{T}, i::TI) where {T <: BlasFloat, TI <: Integer}
    vᵢ = view(V,:,i)
    vᵢ₊₁ = view(V,:,i+1)
    mul!(vᵢ₊₁, A, vᵢ)
    for j = 1:i
        vⱼ = view(V,:,j)
        H[j,i] = dot(vⱼ, vᵢ₊₁)
        axpy!(-H[j,i], vⱼ, vᵢ₊₁)
    end
    β = H[i+1,i] = norm(vᵢ₊₁)
    @. vᵢ₊₁ = vᵢ₊₁ / H[i+1,i]
    return β
end

function arnoldi!(AS::ArnoldiSpace{<:AbstractMatrix{T1}, <:AbstractMatrix{T1}}, A, b::AbstractVector{T2}) where {T1 <: BlasFloat, T2 <: BlasFloat}
    n = size(A, 2)
    V = AS.V
    H = AS.H
    m = AS.m

    n == size(V, 1) || throw(DimensionMismatch())
    n == length(b) || throw(DimensionMismatch())

    arnoldi_init!(A, b, V, H)
    for i = 2:m
        arnoldi_step!(A, V, H, i)
    end
    return AS
end

function arnoldi(A, b::AbstractVector{T}, m::Integer; tol::Real=1e-15) where T <: BlasFloat
    n = size(A, 2)
    V = similar(b, n, m+1)
    H = zeros(T, m+1, m)
    AS = ArnoldiSpace(V, H, m)
    arnoldi!(AS, A, b)
end

### EXPV TOOLS ###

function expv!(x::AbstractVector{T1}, AS::ArnoldiSpace{<:AbstractMatrix{T1}, <:AbstractMatrix{T1}},
    t::T2, b::AbstractVector{T1}) where {T1 <: BlasFloat, T2 <: BlasFloat}

    H = AS.H
    V = AS.V
    m = AS.m

    Hm = view(H, 1:m, 1:m)
    Vm = view(V, :, 1:m)
    lmul!(t, Hm)

    # expH = LinearAlgebra.exp!(Hm)
    expH = exp(Hm)

    cache1 = view(x, 1:m)
    cache2 = view(H, 1, :)
    cache3 = view(H, 2, :)

    mul!(cache1, Vm', b)
    cache2 .= cache1 # just in case cache1 is different from cache2 (e.g., with CUDA)
    mul!(cache3, expH, cache2)
    mul!(x, Vm, cache3)

    return x
end

function expv!(x::AbstractVector{T1}, A, t::T2, b::AbstractVector{T1}; m::Int=min(30, cld(2*length(b), 3))) where {T1 <: BlasFloat, T2 <: BlasFloat}
    AS = arnoldi(A, b, m)
    expv!(x, AS, t, b)
end

function expv(A, t::T1, b::AbstractVector{T2}; m::Int=min(30, cld(2*length(b), 3))) where {T1 <: BlasFloat, T2 <: BlasFloat}
    x = similar(b)
    expv!(x, A, t, b, m=m)
end