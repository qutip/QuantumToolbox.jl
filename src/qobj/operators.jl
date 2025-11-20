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
    rand_unitary([T::Type=ComplexF64,] dimensions, distribution=Val(:haar))

Returns a random unitary [`QuantumObject`](@ref) with element type `T = ComplexF64` (default).

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{ProductDimensions,AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

The `distribution` specifies which of the method used to obtain the unitary matrix:
- `:haar`: Haar random unitary matrix using the algorithm from reference 1
- `:exp`: Uses ``\exp(-i\hat{H})``, where ``\hat{H}`` is a randomly generated Hermitian operator.

# References
1. [F. Mezzadri, How to generate random matrices from the classical compact groups, arXiv:math-ph/0609050 (2007)](https://arxiv.org/abs/math-ph/0609050)

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `rand_unitary(dimensions, Val(distribution))` instead of `rand_unitary(dimensions, distribution)`. Also, put `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
rand_unitary(::Type{T}, dimensions::Int, distribution::Union{Symbol, Val} = Val(:haar)) where {T <: FloatOrComplex} =
    rand_unitary(T, SVector(dimensions), makeVal(distribution))
rand_unitary(::Type{T}, 
    dimensions::Union{ProductDimensions, AbstractVector{Int}, Tuple},
    distribution::Union{Symbol, Val} = Val(:haar),
) where {T <: FloatOrComplex} = rand_unitary(T, dimensions, makeVal(distribution))
function rand_unitary(::Type{T}, dimensions::Union{ProductDimensions, AbstractVector{Int}, Tuple}, ::Val{:haar}) where {T <: FloatOrComplex}
    N = hilbert_dimensions_to_size(dimensions)[1]

    # generate N x N matrix Z of complex standard normal random variates
    Z = randn(T, N, N)

    # find QR decomposition: Z = Q ⋅ R
    Q, R = LinearAlgebra.qr(Z)

    # Create a diagonal matrix Λ by rescaling the diagonal elements of R.
    # Because inv(Λ) ⋅ R has real and strictly positive elements, Q · Λ is therefore Haar distributed.
    Λ = diag(R) # take the diagonal elements of R
    Λ ./= abs.(Λ) # rescaling the elements
    return QuantumObject(to_dense(Q * Diagonal(Λ)); type = Operator(), dims = dimensions)
end
function rand_unitary(::Type{T}, dimensions::Union{ProductDimensions, AbstractVector{Int}, Tuple}, ::Val{:exp}) where {T <: FloatOrComplex}
    N = hilbert_dimensions_to_size(dimensions)[1]

    # generate N x N matrix Z of complex standard normal random variates
    Z = randn(T, N, N)

    # generate Hermitian matrix
    # make it a Qobj first, cause we need Qobj-ver. of exp() for special case in _rand_unitary_exp later
    H = QuantumObject((Z + Z') / 2; type = Operator(), dims = dimensions)
    return to_dense(_rand_unitary_exp(H))
end
rand_unitary(::Type{T}, dimensions::Union{ProductDimensions, AbstractVector{Int}, Tuple}, ::Val{Td}) where {T <: FloatOrComplex, Td} =
    throw(ArgumentError("Invalid distribution: $(Td)"))
rand_unitary(dimensions::Union{Int, Dimensions, AbstractVector{Int}, Tuple}, distribution::Union{Symbol, Val} = Val(:haar)) = rand_unitary(ComplexF64, dimensions, distribution)

# we make H sparse here because the following method currently does not exist : exp(::Matrix{Complex{BigFloat}})
_rand_unitary_exp(H::QuantumObject{Operator, <:AbstractDimensions, Matrix{Complex{BigFloat}}}) = exp(-im * to_sparse(H))
_rand_unitary_exp(H::QuantumObject{Operator, <:AbstractDimensions, <:AbstractArray}) = exp(-im * H)

@doc raw"""
    commutator(A::QuantumObject, B::QuantumObject; anti::Bool=false)

Return the commutator (or `anti`-commutator) of the two [`QuantumObject`](@ref):
- commutator (`anti=false`): ``\hat{A}\hat{B}-\hat{B}\hat{A}``
- anticommutator (`anti=true`): ``\hat{A}\hat{B}+\hat{B}\hat{A}``

Note that `A` and `B` must be [`Operator`](@ref)
"""
commutator(A::QuantumObject{Operator}, B::QuantumObject{Operator}; anti::Bool = false) = A * B - (-1)^anti * B * A

@doc raw"""
    destroy([T::Type=ComplexF64,] N::Int)

Bosonic annihilation operator with Hilbert space cutoff `N` and element type `T = ComplexF64` (default).

This operator acts on a fock state as ``\hat{a} \ket{n} = \sqrt{n} \ket{n-1}``.

# Examples

```jldoctest
julia> a = destroy(20)

Quantum Object:   type=Operator()   dims=[20]   size=(20, 20)   ishermitian=false
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
destroy(::Type{T}, N::Int) where {T <: FloatOrComplex} = QuantumObject(spdiagm(1 => T[sqrt(T(val)) for val in 1:(N - 1)]), Operator(), N)
destroy(N::Int) = destroy(ComplexF64, N)

@doc raw"""
    create([T::Type=ComplexF64,] N::Int)

Bosonic creation operator with Hilbert space cutoff `N` and element type `T = ComplexF64` (default).

This operator acts on a fock state as ``\hat{a}^\dagger \ket{n} = \sqrt{n+1} \ket{n+1}``.

# Examples

```jldoctest
julia> a_d = create(20)

Quantum Object:   type=Operator()   dims=[20]   size=(20, 20)   ishermitian=false
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
create(::Type{T}, N::Int) where {T <: FloatOrComplex} = QuantumObject(spdiagm(-1 => T[sqrt(T(val)) for val in 1:(N - 1)]), Operator(), N)
create(N::Int) = create(ComplexF64, N)

@doc raw"""
    displace(N::Int, α::Number)

Generate a [displacement operator](https://en.wikipedia.org/wiki/Displacement_operator) with element precision same as `α`:

```math
\hat{D}(\alpha)=\exp\left( \alpha \hat{a}^\dagger - \alpha^* \hat{a} \right),
```

where ``\hat{a}`` is the bosonic annihilation operator, and ``\alpha`` is the amount of displacement in optical phase space.
"""
function displace(N::Int, α::T) where {T <: Number}
    a = destroy(_float_type(T), N) # use float type to handle integer α
    return exp(α * a' - conj(α) * a) # the return eltype will be promoted to complex if α is complex
end

@doc raw"""
    squeeze(N::Int, z::Number)

Generate a single-mode [squeeze operator](https://en.wikipedia.org/wiki/Squeeze_operator) with element precision same as `z`:

```math
\hat{S}(z)=\exp\left( \frac{1}{2} (z^* \hat{a}^2 - z(\hat{a}^\dagger)^2) \right),
```

where ``\hat{a}`` is the bosonic annihilation operator.
"""
function squeeze(N::Int, z::T) where {T <: Number}
    a_sq = destroy(_float_type(T), N)^2 # use float type to handle integer z
    return exp((conj(z) * a_sq - z * a_sq') / 2) # the return eltype will be promoted to complex if z is complex
end

@doc raw"""
    num([T::Type=ComplexF64,] N::Int)

Bosonic number operator with Hilbert space cutoff `N` and element type `T = ComplexF64` (default). 

This operator is defined as ``\hat{N}=\hat{a}^\dagger \hat{a}``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
num(::Type{T}, N::Int) where {T <: Number} = QuantumObject(spdiagm(0 => Array{T}(0:(N - 1))), Operator(), N)
num(N::Int) = num(ComplexF64, N)

@doc raw"""
    position([T::Type=ComplexF64,] N::Int)

Position operator with Hilbert space cutoff `N` and element type `T = ComplexF64` (default).

This operator is defined as ``\hat{x}=\frac{1}{\sqrt{2}} (\hat{a}^\dagger + \hat{a})``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
function position(::Type{T}, N::Int) where {T <: FloatOrComplex}
    a = destroy(T, N)
    return (a' + a) / sqrt(T(2))
end
position(N::Int) = position(ComplexF64, N)

@doc raw"""
    momentum([T::Type=ComplexF64,] N::Int)

Momentum operator with Hilbert space cutoff `N` and element type `T = ComplexF64` (default).

This operator is defined as ``\hat{p}= \frac{i}{\sqrt{2}} (\hat{a}^\dagger - \hat{a})``, where ``\hat{a}`` is the bosonic annihilation operator.
"""
function momentum(::Type{T}, N::Int) where {T <: Complex}
    a = destroy(T, N)
    return (a - a') / (im * sqrt(T(2)))
end
momentum(N::Int) = momentum(ComplexF64, N)

@doc raw"""
    phase(N::Int, ϕ0::Real=0)

Single-mode Pegg-Barnett phase operator with Hilbert space cutoff ``N``, the reference phase ``\phi_0``, and element precision same as `ϕ0`.

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
function phase(N::Int, ϕ0::T = 0) where {T <: Real}
    N_list = collect(0:(N - 1))
    ϕ = ϕ0 .+ (2 * _float_type(T)(π) / N) .* N_list # _float_type(T) can deal with T = Int and BigFloat
    states = [exp.(im * ϕ[m] .* N_list) ./ sqrt(T(N)) for m in 1:N] # here promotes element type to complex
    return QuantumObject(sum([ϕ[m] * states[m] * states[m]' for m in 1:N]); type = Operator(), dims = N)
end

@doc raw"""
    jmat([T::Type=ComplexF64,] j::Real, which::Union{Symbol,Val})

Generate higher-order Spin-`j` operators with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

The parameter `which` specifies which of the following operator to return.
- `:x`: ``\hat{S}_x``
- `:y`: ``\hat{S}_y``
- `:z`: ``\hat{S}_z``
- `:+`: ``\hat{S}_+``
- `:-`: ``\hat{S}_-``

Note that if the parameter `which` is not specified, returns a set of Spin-`j` operators: ``(\hat{S}_x, \hat{S}_y, \hat{S}_z)``

# Examples
```jldoctest
julia> jmat(0.5, :x)

Quantum Object:   type=Operator()   dims=[2]   size=(2, 2)   ishermitian=true
2×2 SparseMatrixCSC{ComplexF64, Int64} with 2 stored entries:
     ⋅      0.5+0.0im
 0.5+0.0im      ⋅

julia> jmat(0.5, Val(:-))

Quantum Object:   type=Operator()   dims=[2]   size=(2, 2)   ishermitian=false
2×2 SparseMatrixCSC{ComplexF64, Int64} with 1 stored entry:
     ⋅          ⋅    
 1.0+0.0im      ⋅

julia> jmat(1.5, Val(:z))

Quantum Object:   type=Operator()   dims=[4]   size=(4, 4)   ishermitian=true
4×4 SparseMatrixCSC{ComplexF64, Int64} with 4 stored entries:
 1.5+0.0im      ⋅           ⋅           ⋅    
     ⋅      0.5+0.0im       ⋅           ⋅    
     ⋅          ⋅      -0.5+0.0im       ⋅    
     ⋅          ⋅           ⋅      -1.5+0.0im
```

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `jmat(j, Val(which))` instead of `jmat(j, which)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
jmat(::Type{T}, j::Real, which::Symbol) where {T <: FloatOrComplex} = jmat(T, j, Val(which))
jmat(::Type{T}, j::Real) where {T <: FloatOrComplex} = (jmat(T, j, Val(:x)), jmat(T, j, Val(:y)), jmat(T, j, Val(:z)))
function jmat(::Type{T}, j::Real, ::Val{:x}) where {T <: FloatOrComplex}
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    σ = _jm(T, j)
    return QuantumObject((σ' + σ) / 2, Operator(), Int(J))
end
function jmat(::Type{T}, j::Real, ::Val{:y}) where {T <: FloatOrComplex}
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    σ = _jm(T, j)
    return QuantumObject((σ' - σ) / 2im, Operator(), Int(J))
end
function jmat(::Type{T}, j::Real, ::Val{:z}) where {T <: FloatOrComplex}
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(_jz(T, j), Operator(), Int(J))
end
function jmat(::Type{T}, j::Real, ::Val{:+}) where {T <: FloatOrComplex}
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(adjoint(_jm(T, j)), Operator(), Int(J))
end
function jmat(::Type{T}, j::Real, ::Val{:-}) where {T <: FloatOrComplex}
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    return QuantumObject(_jm(T, j), Operator(), Int(J))
end
jmat(::Type{T}, j::Real, ::Val{T2}) where {T <: FloatOrComplex, T2} = throw(ArgumentError("Invalid spin operator: $(T2)"))
jmat(j::Real, which::Union{Symbol, Val}) = jmat(ComplexF64, j, which)
jmat(j::Real) = jmat(ComplexF64, j)

function _jm(::Type{T}, j::Real) where {T <: FloatOrComplex}
    m = T.(j:(-1):(-j))
    data = @. sqrt(T(j * (j + 1)) - m * (m - 1))
    return spdiagm(-1 => data[1:(end - 1)])
end
function _jz(::Type{T}, j::Real) where {T <: FloatOrComplex}
    data = @. T(j) - (0:Int(2 * j))
    return spdiagm(0 => data)
end

@doc raw"""
    spin_Jx([T::Type=ComplexF64,] j::Real)

``\hat{S}_x`` operator for Spin-`j` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

See also [`jmat`](@ref).
"""
spin_Jx(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j, Val(:x))
spin_Jx(j::Real) = spin_Jx(ComplexF64, j)

@doc raw"""
    spin_Jy([T::Type=ComplexF64,] j::Real)

``\hat{S}_y`` operator for Spin-`j` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

See also [`jmat`](@ref).
"""
spin_Jy(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j, Val(:y))
spin_Jy(j::Real) = spin_Jy(ComplexF64, j)

@doc raw"""
    spin_Jz([T::Type=ComplexF64,] j::Real)

``\hat{S}_z`` operator for Spin-`j` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

See also [`jmat`](@ref).
"""
spin_Jz(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j, Val(:z))
spin_Jz(j::Real) = spin_Jz(ComplexF64, j)

@doc raw"""
    spin_Jm([T::Type=ComplexF64,] j::Real)

``\hat{S}_-`` operator for Spin-`j` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

See also [`jmat`](@ref).
"""
spin_Jm(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j, Val(:-))
spin_Jm(j::Real) = spin_Jm(ComplexF64, j)

@doc raw"""
    spin_Jp([T::Type=ComplexF64,] j::Real)

``\hat{S}_+`` operator for Spin-`j` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

See also [`jmat`](@ref).
"""
spin_Jp(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j, Val(:+))
spin_Jp(j::Real) = spin_Jp(ComplexF64, j)

@doc raw"""
    spin_J_set([T::Type=ComplexF64,] j::Real)

A set of Spin-`j` operators ``(\hat{S}_x, \hat{S}_y, \hat{S}_z)`` with element type `T = ComplexF64` (default). `j` is the spin quantum number and can be a non-negative integer or half-integer.

Note that this functions is same as `jmat(j)`. See also [`jmat`](@ref).
"""
spin_J_set(::Type{T}, j::Real) where {T <: FloatOrComplex} = jmat(T, j)
spin_J_set(j::Real) = spin_J_set(ComplexF64, j)

@doc raw"""
    sigmap([T::Type=ComplexF64])

Pauli ladder operator ``\hat{\sigma}_+ = (\hat{\sigma}_x + i \hat{\sigma}_y) / 2`` with element type `T = ComplexF64` (default).

See also [`jmat`](@ref).
"""
sigmap(::Type{T}) where {T <: FloatOrComplex} = jmat(T, 0.5, Val(:+))
sigmap() = sigmap(ComplexF64)

@doc raw"""
    sigmam([T::Type=ComplexF64])

Pauli ladder operator ``\hat{\sigma}_- = (\hat{\sigma}_x - i \hat{\sigma}_y) / 2`` with element type `T = ComplexF64` (default).

See also [`jmat`](@ref).
"""
sigmam(::Type{T}) where {T <: FloatOrComplex} = jmat(T, 0.5, Val(:-))
sigmam() = sigmam(ComplexF64)

@doc raw"""
    sigmax([T::Type=ComplexF64])

Pauli operator ``\hat{\sigma}_x = \hat{\sigma}_- + \hat{\sigma}_+`` with element type `T = ComplexF64` (default).

See also [`jmat`](@ref).
"""
sigmax(::Type{T}) where {T <: FloatOrComplex} = rmul!(jmat(T, 0.5, Val(:x)), 2)
sigmax() = sigmax(ComplexF64)

@doc raw"""
    sigmay([T::Type=ComplexF64])

Pauli operator ``\hat{\sigma}_y = i \left( \hat{\sigma}_- - \hat{\sigma}_+ \right)`` with element type `T = ComplexF64` (default).

See also [`jmat`](@ref).
"""
sigmay(::Type{T}) where {T <: FloatOrComplex} = rmul!(jmat(T, 0.5, Val(:y)), 2)
sigmay() = sigmay(ComplexF64)

@doc raw"""
    sigmaz([T::Type=ComplexF64])

Pauli operator ``\hat{\sigma}_z = \left[ \hat{\sigma}_+ , \hat{\sigma}_- \right]`` with element type `T = ComplexF64` (default).

See also [`jmat`](@ref).
"""
sigmaz(::Type{T}) where {T <: FloatOrComplex} = rmul!(jmat(T, 0.5, Val(:z)), 2)
sigmaz() = sigmaz(ComplexF64)

@doc raw"""
    eye([T::Type=ComplexF64,] N::Int; type=Operator, dims=nothing)
    qeye([T::Type=ComplexF64,] N::Int; type=Operator, dims=nothing)

Identity operator ``\hat{\mathbb{1}}`` with size `N` and element type `T = ComplexF64` (default).

It is also possible to specify the list of Hilbert dimensions `dims` if different subsystems are present.

Note that `type` can only be either [`Operator`](@ref) or [`SuperOperator`](@ref)

!!! note
    `qeye` is a synonym of `eye`.
"""
function eye(::Type{T}, N::Int; type = Operator(), dims = nothing) where {T <: Number}
    if dims isa Nothing
        dims = isa(type, Operator) ? N : isqrt(N)
    end
    return QuantumObject(Diagonal(ones(T, N)); type, dims)
end
eye(N::Int; type = Operator(), dims = nothing) = eye(ComplexF64, N; type, dims)

@doc raw"""
    fdestroy([T::Type=ComplexF64,] N::Union{Int,Val}, j::Int)

Construct a fermionic destruction operator [with element type `T = ComplexF64` (default)] acting on the `j`-th site, where the fock space has totally `N`-sites:

Here, we use the [Jordan-Wigner transformation](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation), namely
```math
\hat{d}_j = \hat{\sigma}_z^{\otimes j-1} \otimes \hat{\sigma}_{+} \otimes \hat{\mathbb{1}}^{\otimes N-j}
```

The site index `j` should satisfy: `1 ≤ j ≤ N`.

Note that we put ``\hat{\sigma}_{+} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}`` here because we consider ``|0\rangle = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`` to be ground (vacant) state, and ``|1\rangle = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`` to be excited (occupied) state.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `fdestroy(Val(N), j)` instead of `fdestroy(N, j)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
fdestroy(::Type{T}, N::Union{Int, Val}, j::Int) where {T <: FloatOrComplex} = _Jordan_Wigner(T, N, j, sigmap(T))
fdestroy(N::Union{Int, Val}, j::Int) = fdestroy(ComplexF64, N, j)

@doc raw"""
    fcreate([T::Type=ComplexF64,] N::Union{Int,Val}, j::Int)

Construct a fermionic creation operator [with element type `T = ComplexF64` (default)] acting on the `j`-th site, where the fock space has totally `N`-sites:

Here, we use the [Jordan-Wigner transformation](https://en.wikipedia.org/wiki/Jordan%E2%80%93Wigner_transformation), namely
```math
\hat{d}^\dagger_j = \hat{\sigma}_z^{\otimes j-1} \otimes \hat{\sigma}_{-} \otimes \hat{\mathbb{1}}^{\otimes N-j}
```

The site index `j` should satisfy: `1 ≤ j ≤ N`.

Note that we put ``\hat{\sigma}_{-} = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}`` here because we consider ``|0\rangle = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`` to be ground (vacant) state, and ``|1\rangle = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`` to be excited (occupied) state.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `fcreate(Val(N), j)` instead of `fcreate(N, j)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
fcreate(::Type{T}, N::Union{Int, Val}, j::Int) where {T <: FloatOrComplex} = _Jordan_Wigner(T, N, j, sigmam(T))
fcreate(N::Union{Int, Val}, j::Int) = fcreate(ComplexF64, N, j)

_Jordan_Wigner(::Type{T}, N::Int, j::Int, op::QuantumObject{Operator}) where {T <: FloatOrComplex} = _Jordan_Wigner(T, Val(N), j, op)

function _Jordan_Wigner(::Type{T}, ::Val{N}, j::Int, op::QuantumObject{Operator}) where {T <: FloatOrComplex, N}
    (N < 1) && throw(ArgumentError("The total number of sites (N) cannot be less than 1"))
    (1 <= j <= N) || throw(ArgumentError("The site index (j) should satisfy: 1 ≤ j ≤ N"))

    σz = sigmaz(T).data
    Z_tensor = kron(one(T), one(T), fill(σz, j - 1)...)

    S = 2^(N - j)
    I_tensor = sparse(one(T) * LinearAlgebra.I, S, S)

    return QuantumObject(kron(Z_tensor, op.data, I_tensor); type = Operator(), dims = ntuple(i -> 2, Val(N)))
end

@doc raw"""
    projection([T::Type=ComplexF64,] N::Int, i::Int, j::Int)

Generates the projection operator ``\hat{O} = |i \rangle\langle j|`` with Hilbert space dimension `N` and element type `T = ComplexF64` (default).
"""
function projection(::Type{T}, N::Int, i::Int, j::Int) where {T <: Number}
    (0 <= i < N) || throw(ArgumentError("Invalid argument i, must satisfy: 0 ≤ i ≤ N-1"))
    (0 <= j < N) || throw(ArgumentError("Invalid argument j, must satisfy: 0 ≤ j ≤ N-1"))

    return QuantumObject(sparse([i + 1], [j + 1], [one(T)], N, N), type = Operator(), dims = N)
end
projection(N::Int, i::Int, j::Int) = projection(ComplexF64, N, i, j)

@doc raw"""
    tunneling([T::Type=ComplexF64,] N::Int, m::Int=1; sparse::Union{Bool,Val{<:Bool}}=Val(false))

Generate a tunneling operator with element type `T = ComplexF64` (default), defined as:

```math
\sum_{n=0}^{N-m} | n \rangle\langle n+m | + | n+m \rangle\langle n |,
```

where ``N`` is the number of basis states in the Hilbert space, and ``m`` is the number of excitations in tunneling event.

If `sparse=true`, the operator is returned as a sparse matrix, otherwise a dense matrix is returned.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `tunneling(N, m, Val(sparse))` instead of `tunneling(N, m, sparse)`. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function tunneling(::Type{T}, N::Int, m::Int = 1; sparse::Union{Bool, Val} = Val(false)) where {T <: Number}
    (m < 1) && throw(ArgumentError("The number of excitations (m) cannot be less than 1"))

    data = ones(T, N - m)
    if getVal(sparse)
        return QuantumObject(spdiagm(m => data, -m => data); type = Operator(), dims = N)
    else
        return QuantumObject(diagm(m => data, -m => data); type = Operator(), dims = N)
    end
end
tunneling(N::Int, m::Int = 1; sparse::Union{Bool, Val} = Val(false)) = tunneling(ComplexF64, N, m; sparse)

@doc raw"""
    qft([T::Type=ComplexF64,] dimensions)

Generates a discrete Fourier transform matrix ``\hat{F}_N`` for [Quantum Fourier Transform (QFT)](https://en.wikipedia.org/wiki/Quantum_Fourier_transform) with given argument `dimensions` and element type `T = ComplexF64` (default).

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Union{ProductDimensions,AbstractVector{Int},Tuple}`: list of dimensions representing the each number of basis in the subsystems.

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

!!! warning "Beware of type-stability!"
    It is highly recommended to use `qft(dimensions)` with `dimensions` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
qft(::Type{T}, dimensions::Int) where {T <: Complex} = QuantumObject(_qft_op(T, dimensions), Operator(), dimensions)
qft(::Type{T}, dimensions::Union{ProductDimensions, AbstractVector{Int}, Tuple}) where {T <: Complex} =
    QuantumObject(_qft_op(T, hilbert_dimensions_to_size(dimensions)[1]), Operator(), dimensions)
qft(dimensions::Union{Int, Dimensions, AbstractVector{Int}, Tuple}) = qft(ComplexF64, dimensions)

function _qft_op(::Type{T}, N::Int) where {T <: Complex}
    N_T = T(N)
    ω = exp(2 * im * T(π) / N_T)
    arr = 0:(N - 1)
    L, M = meshgrid(arr, arr)
    return ω .^ (L .* M) / sqrt(N_T)
end
