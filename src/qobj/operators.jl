#=
Functions for generating (common) quantum operators.
=#

export rand_unitary
export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set
export sigmam, sigmap, sigmax, sigmay, sigmaz
export destroy, create, eye, projection
export displace, squeeze, num, phase
export fdestroy, fcreate
export commutator
export tunneling
export qft

@doc raw"""
    rand_unitary(dimensions, distribution=:haar)

Returns a random unitary [`QuantumObject`](@ref).

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

The `distribution` specifies which of the method used to obtain the unitary matrix:
- `:haar`: Haar random unitary matrix using the algorithm from reference 1
- `:exp`: Uses ``\exp(-iH)``, where ``H`` is a randomly generated Hermitian operator.

# References
1. [F. Mezzadri, How to generate random matrices from the classical compact groups, arXiv:math-ph/0609050 (2007)](https://arxiv.org/abs/math-ph/0609050)
"""
rand_unitary(dimensions::Int, distribution::Symbol = :haar) = rand_unitary(SVector(dimensions), Val(distribution))
rand_unitary(dimensions::Union{AbstractVector{Int},Tuple}, distribution::Symbol = :haar) =
    rand_unitary(dimensions, Val(distribution))
function rand_unitary(dimensions::Union{AbstractVector{Int},Tuple}, ::Val{:haar})
    N = prod(dimensions)

    # generate N x N matrix Z of complex standard normal random variates
    Z = randn(ComplexF64, N, N)

    # find QR decomposition: Z = Q ⋅ R
    Q, R = LinearAlgebra.qr(Z)

    # Create a diagonal matrix Λ by rescaling the diagonal elements of R.
    # Because inv(Λ) ⋅ R has real and strictly positive elements, Q · Λ is therefore Haar distributed.
    Λ = diag(R) # take the diagonal elements of R
    Λ ./= abs.(Λ) # rescaling the elements
    return QuantumObject(dense_to_sparse(Q * Diagonal(Λ)); type = Operator, dims = dimensions)
end
function rand_unitary(dimensions::Union{AbstractVector{Int},Tuple}, ::Val{:exp})
    N = prod(dimensions)

    # generate N x N matrix Z of complex standard normal random variates
    Z = randn(ComplexF64, N, N)

    # generate Hermitian matrix
    H = QuantumObject((Z + Z') / 2; type = Operator, dims = dimensions)

    return exp(-1.0im * H)
end
rand_unitary(dimensions::Union{AbstractVector{Int},Tuple}, ::Val{T}) where {T} =
    throw(ArgumentError("Invalid distribution: $(T)"))

@doc raw"""
    commutator(A::QuantumObject, B::QuantumObject; anti::Bool=false)

Return the commutator (or `anti`-commutator) of the two [`QuantumObject`](@ref):
- commutator (`anti=false`): ``AB-BA``
- anticommutator (`anti=true`): ``AB+BA``

Note that `A` and `B` must be [`Operator`](@ref)
"""
commutator(
    A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject};
    anti::Bool = false,
) where {T1,T2} = A * B - (-1)^anti * B * A

@doc raw"""
    destroy(N::Int)

Bosonic annihilation operator with Hilbert space cutoff `N`. 

This operator acts on a fock state as ``\hat{a} \ket{n} = \sqrt{n} \ket{n-1}``.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⎡⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⎦

julia> fock(20, 3)' * a * fock(20, 4)
2.0 + 0.0im
```
"""
destroy(N::Int) = QuantumObject(spdiagm(1 => Array{ComplexF64}(sqrt.(1:N-1))), Operator, N)

@doc raw"""
    create(N::Int)

Bosonic creation operator with Hilbert space cutoff `N`.

This operator acts on a fock state as ``\hat{a}^\dagger \ket{n} = \sqrt{n+1} \ket{n+1}``.

# Examples

```
julia> a_d = create(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⎡⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠢⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠢⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠢⡀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⎦

julia> fock(20, 4)' * a_d * fock(20, 3)
2.0 + 0.0im
```
"""
create(N::Int) = QuantumObject(spdiagm(-1 => Array{ComplexF64}(sqrt.(1:N-1))), Operator, N)

@doc raw"""
    displace(N::Int, α::Number)

Generate a [displacement operator](https://en.wikipedia.org/wiki/Displacement_operator):

```math
\hat{D}(\alpha)=\exp\left( \alpha \hat{a}^\dagger - \alpha^* \hat{a} \right),
```

where ``\hat{a}`` is the bosonic annihilation operator, and ``\alpha`` is the amount of displacement in optical phase space.
"""
function displace(N::Int, α::T) where {T<:Number}
    a = destroy(N)
    return exp(α * a' - α' * a)
end

@doc raw"""
    squeeze(N::Int, z::Number)

Generate a single-mode [squeeze operator](https://en.wikipedia.org/wiki/Squeeze_operator):

```math
\hat{S}(z)=\exp\left( \frac{1}{2} (z^* \hat{a}^2 - z(\hat{a}^\dagger)^2) \right),
```

where ``\hat{a}`` is the bosonic annihilation operator.
"""
function squeeze(N::Int, z::T) where {T<:Number}
    a_sq = destroy(N)^2
    return exp((z' * a_sq - z * a_sq') / 2)
end

@doc raw"""
    num(N::Int)

Bosonic number operator with Hilbert space cutoff `N`. 

This operator is defined as ``\hat{N}=\hat{a}^\dagger \hat{a}``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
num(N::Int) = QuantumObject(spdiagm(0 => Array{ComplexF64}(0:N-1)), Operator, N)

@doc raw"""
    position(N::Int)

Position operator with Hilbert space cutoff `N`. 

This operator is defined as ``\hat{x}=\frac{1}{\sqrt{2}} (\hat{a}^\dagger + \hat{a})``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
function position(N::Int)
    a = destroy(N)
    return (a' + a) / sqrt(2)
end

@doc raw"""
    momentum(N::Int)

Momentum operator with Hilbert space cutoff `N`. 

This operator is defined as ``\hat{p}= \frac{i}{\sqrt{2}} (\hat{a}^\dagger - \hat{a})``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
function momentum(N::Int)
    a = destroy(N)
    return (1.0im * sqrt(0.5)) * (a' - a)
end

@doc raw"""
    phase(N::Int, ϕ0::Real=0)

Single-mode Pegg-Barnett phase operator with Hilbert space cutoff ``N`` and the reference phase ``\phi_0``.

This operator is defined as

```math
\hat{\phi} = \sum_{m=0}^{N-1} \phi_m |\phi_m\rangle \langle\phi_m|,
```

where

```math
\phi_m = \phi_0 + \frac{2m\pi}{N},
```

and

```math
|\phi_m\rangle = \frac{1}{\sqrt{N}} \sum_{n=0}^{N-1} \exp(i n \phi_m) |n\rangle.
```

# Reference
- [Michael Martin Nieto, QUANTUM PHASE AND QUANTUM PHASE OPERATORS: Some Physics and Some History, arXiv:hep-th/9304036](https://arxiv.org/abs/hep-th/9304036), Equation (30-32).
"""
function phase(N::Int, ϕ0::Real = 0)
    N_list = collect(0:(N-1))
    ϕ = ϕ0 .+ (2 * π / N) .* N_list
    states = [exp.((1.0im * ϕ[m]) .* N_list) ./ sqrt(N) for m in 1:N]
    return QuantumObject(sum([ϕ[m] * states[m] * states[m]' for m in 1:N]); type = Operator, dims = N)
end

@doc raw"""
    jmat(j::Real, which::Symbol)

Generate higher-order Spin-`j` operators, where `j` is the spin quantum number and can be a non-negative integer or half-integer

The parameter `which` specifies which of the following operator to return.
- `:x`: ``S_x``
- `:y`: ``S_y``
- `:z`: ``S_z``
- `:+`: ``S_+``
- `:-`: ``S_-``

Note that if the parameter `which` is not specified, returns a set of Spin-`j` operators: ``(S_x, S_y, S_z)``

# Examples
```
julia> jmat(0.5, :x)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 SparseMatrixCSC{ComplexF64, Int64} with 2 stored entries:
     ⋅      0.5+0.0im
 0.5+0.0im      ⋅

julia> jmat(0.5, :-)
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false
2×2 SparseMatrixCSC{ComplexF64, Int64} with 1 stored entry:
     ⋅          ⋅    
 1.0+0.0im      ⋅
```
"""
jmat(j::Real, which::Symbol) = jmat(j, Val(which))
jmat(j::Real) = (jmat(j, Val(:x)), jmat(j, Val(:y)), jmat(j, Val(:z)))
function jmat(j::Real, ::Val{:x})
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    σ = _jm(j)
    return QuantumObject((σ' + σ) / 2, Operator, Int(J))
end
function jmat(j::Real, ::Val{:y})
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    σ = _jm(j)
    return QuantumObject((σ' - σ) / 2im, Operator, Int(J))
end
function jmat(j::Real, ::Val{:z})
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(_jz(j), Operator, Int(J))
end
function jmat(j::Real, ::Val{:+})
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(adjoint(_jm(j)), Operator, Int(J))
end
function jmat(j::Real, ::Val{:-})
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(_jm(j), Operator, Int(J))
end
jmat(j::Real, ::Val{T}) where {T} = throw(ArgumentError("Invalid spin operator: $(T)"))

function _jm(j::Real)
    m = j:(-1):-j
    data = sqrt.(j * (j + 1) .- m .* (m .- 1))[1:end-1]
    return spdiagm(-1 => Array{ComplexF64}(data))
end
_jz(j::Real) = spdiagm(0 => Array{ComplexF64}(j .- (0:Int(2 * j))))

@doc raw"""
    spin_Jx(j::Real)

``S_x`` operator for Spin-`j`, where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
spin_Jx(j::Real) = jmat(j, Val(:x))

@doc raw"""
    spin_Jy(j::Real)

``S_y`` operator for Spin-`j`, where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
spin_Jy(j::Real) = jmat(j, Val(:y))

@doc raw"""
    spin_Jz(j::Real)

``S_z`` operator for Spin-`j`, where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
spin_Jz(j::Real) = jmat(j, Val(:z))

@doc raw"""
    spin_Jm(j::Real)

``S_-`` operator for Spin-`j`, where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
spin_Jm(j::Real) = jmat(j, Val(:-))

@doc raw"""
    spin_Jp(j::Real)

``S_+`` operator for Spin-`j`, where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
spin_Jp(j::Real) = jmat(j, Val(:+))

@doc raw"""
    spin_J_set(j::Real)

A set of Spin-`j` operators ``(S_x, S_y, S_z)``, where `j` is the spin quantum number and can be a non-negative integer or half-integer

Note that this functions is same as `jmat(j)`. See also [`jmat`](@ref).
"""
spin_J_set(j::Real) = jmat(j)

@doc raw"""
    sigmap()

Pauli ladder operator ``\hat{\sigma}_+ = \hat{\sigma}_x + i \hat{\sigma}_y``.

See also [`jmat`](@ref).
"""
sigmap() = jmat(0.5, Val(:+))

@doc raw"""
    sigmam()

Pauli ladder operator ``\hat{\sigma}_- = \hat{\sigma}_x - i \hat{\sigma}_y``.

See also [`jmat`](@ref).
"""
sigmam() = jmat(0.5, Val(:-))

@doc raw"""
    sigmax()

Pauli operator ``\hat{\sigma}_x = \hat{\sigma}_- + \hat{\sigma}_+``.

See also [`jmat`](@ref).
"""
sigmax() = rmul!(jmat(0.5, Val(:x)), 2)

@doc raw"""
    sigmay()

Pauli operator ``\hat{\sigma}_y = i \left( \hat{\sigma}_- - \hat{\sigma}_+ \right)``.

See also [`jmat`](@ref).
"""
sigmay() = rmul!(jmat(0.5, Val(:y)), 2)

@doc raw"""
    sigmaz()

Pauli operator ``\hat{\sigma}_z = \comm{\hat{\sigma}_+}{\hat{\sigma}_-}``.

See also [`jmat`](@ref).
"""
sigmaz() = rmul!(jmat(0.5, Val(:z)), 2)

@doc raw"""
    eye(N::Int; type=Operator, dims=nothing)

Identity operator ``\hat{\mathbb{1}}`` with size `N`.

It is also possible to specify the list of Hilbert dimensions `dims` if different subsystems are present.

Note that `type` can only be either [`Operator`](@ref) or [`SuperOperator`](@ref)
"""
eye(
    N::Int;
    type::ObjType = Operator,
    dims = nothing,
) where {ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Diagonal(ones(ComplexF64, N)); type = type, dims = dims)

@doc raw"""
    fdestroy(N::Union{Int,Val}, j::Int)

Construct a fermionic destruction operator acting on the `j`-th site, where the fock space has totally `N`-sites:

Here, we use the [Jordan-Wigner transformation](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation), namely
```math
d_j = \sigma_z^{\otimes j} \otimes \sigma_{-} \otimes I^{\otimes N-j-1}
```

Note that the site index `j` should satisfy: `0 ≤ j ≤ N - 1`.

> [!IMPORTANT]
> If you want to keep type stability, it is recommended to use `fdestroy(Val(N), j)` instead of `fdestroy(N, j)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
fdestroy(N::Union{Int,Val}, j::Int) = _Jordan_Wigner(N, j, sigmam())

@doc raw"""
    fcreate(N::Union{Int,Val}, j::Int)

Construct a fermionic creation operator acting on the `j`-th site, where the fock space has totally `N`-sites:

Here, we use the [Jordan-Wigner transformation](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation), namely
```math
d_j^\dagger = \sigma_z^{\otimes j} \otimes \sigma_{+} \otimes I^{\otimes N-j-1}
```

Note that the site index `j` should satisfy: `0 ≤ j ≤ N - 1`.

> [!IMPORTANT]
> If you want to keep type stability, it is recommended to use `fcreate(Val(N), j)` instead of `fcreate(N, j)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
fcreate(N::Union{Int,Val}, j::Int) = _Jordan_Wigner(N, j, sigmap())

_Jordan_Wigner(N::Int, j::Int, op::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} =
    _Jordan_Wigner(Val(N), j, op)

function _Jordan_Wigner(::Val{N}, j::Int, op::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {N,T}
    (N < 1) && throw(ArgumentError("The total number of sites (N) cannot be less than 1"))
    ((j >= N) || (j < 0)) && throw(ArgumentError("The site index (j) should satisfy: 0 ≤ j ≤ N - 1"))

    σz = sigmaz().data
    Z_tensor = kron(1, 1, fill(σz, j)...)

    S = 2^(N - j - 1)
    I_tensor = sparse((1.0 + 0.0im) * LinearAlgebra.I, S, S)

    return QuantumObject(kron(Z_tensor, op.data, I_tensor); type = Operator, dims = ntuple(i -> 2, Val(N)))
end

@doc raw"""
    projection(N::Int, i::Int, j::Int)

Generates the projection operator ``\hat{O} = \dyad{i}{j}`` with Hilbert space dimension `N`.
"""
projection(N::Int, i::Int, j::Int) = QuantumObject(sparse([i + 1], [j + 1], [1.0 + 0.0im], N, N))

@doc raw"""
    tunneling(N::Int, m::Int=1; sparse::Union{Bool,Val{<:Bool}}=Val(false))

Generate a tunneling operator defined as:

```math
\sum_{n=0}^{N-m} | n \rangle\langle n+m | + | n+m \rangle\langle n |,
```

where ``N`` is the number of basis states in the Hilbert space, and ``m`` is the number of excitations in tunneling event.

If `sparse=true`, the operator is returned as a sparse matrix, otherwise a dense matrix is returned.

> [!IMPORTANT]
> If you want to keep type stability, it is recommended to use `tunneling(N, m, Val(sparse))` instead of `tunneling(N, m, sparse)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) for more details.
"""
function tunneling(N::Int, m::Int = 1; sparse::Union{Bool,Val} = Val(false))
    (m < 1) && throw(ArgumentError("The number of excitations (m) cannot be less than 1"))

    data = ones(ComplexF64, N - m)
    if getVal(makeVal(sparse))
        return QuantumObject(spdiagm(m => data, -m => data); type = Operator, dims = N)
    else
        return QuantumObject(diagm(m => data, -m => data); type = Operator, dims = N)
    end
end

@doc raw"""
    qft(dimensions)

Generates a discrete Fourier transform matrix ``\hat{F}_N`` for [Quantum Fourier Transform (QFT)](https://en.wikipedia.org/wiki/Quantum_Fourier_transform) with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

``N`` represents the total dimension, and therefore the matrix is defined as

```math
\hat{F}_N = \frac{1}{\sqrt{N}}\begin{bmatrix}
1 & 1 & 1 & 1 & \cdots & 1\\
1 & \omega & \omega^2 & \omega^3 & \cdots & \omega^{N-1}\\
1 & \omega^2 & \omega^4 & \omega^6 & \cdots & \omega^{2(N-1)}\\
1 & \omega^3 & \omega^6 & \omega^9 & \cdots & \omega^{3(N-1)}\\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots\\
1 & \omega^{N-1} & \omega^{2(N-1)} & \omega^{3(N-1)} & \cdots & \omega^{(N-1)(N-1)}
\end{bmatrix},
```

where ``\omega = \exp(\frac{2 \pi i}{N})``.
"""
qft(dimensions::Int) = QuantumObject(_qft_op(dimensions), Operator, dimensions)
qft(dimensions::Union{AbstractVector{T},Tuple}) where {T} =
    QuantumObject(_qft_op(prod(dimensions)), Operator, dimensions)
function _qft_op(N::Int)
    ω = exp(2.0im * π / N)
    arr = 0:(N-1)
    L, M = meshgrid(arr, arr)
    return ω .^ (L .* M) / sqrt(N)
end
