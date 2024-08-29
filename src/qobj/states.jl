#=
Functions for generating (common) quantum states.
=#

export zero_ket, fock, basis, coherent, rand_ket
export fock_dm, coherent_dm, thermal_dm, maximally_mixed_dm, rand_dm
export spin_state, spin_coherent
export bell_state, singlet_state, triplet_states, w_state, ghz_state

@doc raw"""
    zero_ket(dimensions)

Returns a zero [`Ket`](@ref) vector with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.
"""
zero_ket(dimensions::Int) = QuantumObject(zeros(ComplexF64, dimensions), Ket, [dimensions])
zero_ket(dimensions::Vector{Int}) = QuantumObject(zeros(ComplexF64, prod(dimensions)), Ket, dimensions)

@doc raw"""
    fock(N::Int, pos::Int=0; dims::Vector{Int}=[N], sparse::Bool=false)

Generates a fock state ``\ket{\psi}`` of dimension `N`. 

It is also possible to specify the list of dimensions `dims` if different subsystems are present.
"""
function fock(N::Int, pos::Int = 0; dims::Vector{Int} = [N], sparse::Bool = false)
    if sparse
        array = sparsevec([pos + 1], [1.0 + 0im], N)
    else
        array = zeros(ComplexF64, N)
        array[pos+1] = 1
    end
    return QuantumObject(array; type = Ket, dims = dims)
end

@doc raw"""
    basis(N::Int, pos::Int = 0; dims::Vector{Int}=[N])

Generates a fock state like [`fock`](@ref).

It is also possible to specify the list of dimensions `dims` if different subsystems are present.
"""
basis(N::Int, pos::Int = 0; dims::Vector{Int} = [N]) = fock(N, pos, dims = dims)

@doc raw"""
    coherent(N::Int, α::Number)

Generates a [coherent state](https://en.wikipedia.org/wiki/Coherent_state) ``|\alpha\rangle``, which is defined as an eigenvector of the bosonic annihilation operator ``\hat{a} |\alpha\rangle = \alpha |\alpha\rangle``.

This state is constructed via the displacement operator [`displace`](@ref) and zero-fock state [`fock`](@ref): ``|\alpha\rangle = \hat{D}(\alpha) |0\rangle``
"""
coherent(N::Int, α::T) where {T<:Number} = displace(N, α) * fock(N, 0)

@doc raw"""
    rand_ket(dimensions)

Generate a random normalized [`Ket`](@ref) vector with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.
"""
rand_ket(dimensions::Int) = rand_ket([dimensions])
function rand_ket(dimensions::Vector{Int})
    N = prod(dimensions)
    ψ = rand(ComplexF64, N) .- (0.5 + 0.5im)
    return QuantumObject(normalize!(ψ); type = Ket, dims = dimensions)
end

@doc raw"""
    fock_dm(N::Int, pos::Int=0; dims::Vector{Int}=[N], sparse::Bool=false)

Density matrix representation of a Fock state.

Constructed via outer product of [`fock`](@ref).
"""
function fock_dm(N::Int, pos::Int = 0; dims::Vector{Int} = [N], sparse::Bool = false)
    ψ = fock(N, pos; dims = dims, sparse = sparse)
    return ψ * ψ'
end

@doc raw"""
    coherent_dm(N::Int, α::Number)

Density matrix representation of a [coherent state](https://en.wikipedia.org/wiki/Coherent_state).

Constructed via outer product of [`coherent`](@ref).
"""
function coherent_dm(N::Int, α::T) where {T<:Number}
    ψ = coherent(N, α)
    return ψ * ψ'
end

@doc raw"""
    thermal_dm(N::Int, n::Real; sparse::Bool=false)

Density matrix for a thermal state (generating thermal state probabilities) with the following arguments:
- `N::Int`: Number of basis states in the Hilbert space
- `n::Real`: Expectation value for number of particles in the thermal state.
"""
function thermal_dm(N::Int, n::Real; sparse::Bool = false)
    β = log(1.0 / n + 1.0)
    N_list = Array{Float64}(0:N-1)
    data = exp.(-β .* N_list)
    if sparse
        return QuantumObject(spdiagm(0 => data ./ sum(data)), Operator, [N])
    else
        return QuantumObject(diagm(0 => data ./ sum(data)), Operator, [N])
    end
end

@doc raw"""
    maximally_mixed_dm(dimensions)

Returns the maximally mixed density matrix with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.
"""
maximally_mixed_dm(dimensions::Int) = QuantumObject(I(dimensions) / complex(dimensions), Operator, [dimensions])
function maximally_mixed_dm(dimensions::Vector{Int})
    N = prod(dimensions)
    return QuantumObject(I(N) / complex(N), Operator, dimensions)
end

@doc raw"""
    rand_dm(dimensions; rank::Int=prod(dimensions))

Generate a random density matrix from Ginibre ensemble with given argument `dimensions` and `rank`, ensuring that it is positive semi-definite and trace equals to `1`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.

The default keyword argument `rank = prod(dimensions)` (full rank).

# References
- [J. Ginibre, Statistical ensembles of complex, quaternion, and real matrices, Journal of Mathematical Physics 6.3 (1965): 440-449](https://doi.org/10.1063/1.1704292)
- [K. Życzkowski, et al., Generating random density matrices, Journal of Mathematical Physics 52, 062201 (2011)](http://dx.doi.org/10.1063/1.3595693)
"""
rand_dm(dimensions::Int; rank::Int = prod(dimensions)) = rand_dm([dimensions], rank = rank)
function rand_dm(dimensions::Vector{Int}; rank::Int = prod(dimensions))
    N = prod(dimensions)
    (rank < 1) && throw(DomainError(rank, "The argument rank must be larger than 1."))
    (rank > N) && throw(DomainError(rank, "The argument rank cannot exceed dimensions."))

    X = _Ginibre_ensemble(N, rank)
    ρ = X * X'
    ρ /= tr(ρ)
    return QuantumObject(ρ; type = Operator, dims = dimensions)
end

@doc raw"""
    spin_state(j::Real, m::Real)

Generate the spin state: ``|j, m\rangle``

The eigenstate of the Spin-`j` ``S_z`` operator with eigenvalue `m`, where where `j` is the spin quantum number and can be a non-negative integer or half-integer

See also [`jmat`](@ref).
"""
function spin_state(j::Real, m::Real)
    J = 2 * j + 1
    ((floor(J) != J) || (j < 0)) &&
        throw(ArgumentError("The spin quantum number (j) must be a non-negative integer or half-integer."))

    Δ = j - m
    ((floor(Δ) != Δ) || (Δ < 0)) &&
        throw(ArgumentError("Invalid eigenvalue m: (j - m) must be a non-negative integer."))
    (m < (-j)) && throw(ArgumentError("Invalid eigenvalue m, must satisfy: -j ≤ m ≤ j"))

    return basis(Int(J), Int(Δ))
end

@doc raw"""
    spin_coherent(j::Real, θ::Real, ϕ::Real)

Generate the coherent spin state (rotation of the ``|j, j\rangle`` state), namely

```math
|\theta, \phi \rangle = R(\theta, \phi) |j, j\rangle
```

where the rotation operator is defined as

```math
R(\theta, \phi) = \exp \left( \frac{\theta}{2} (S_- e^{i\phi} - S_+ e^{-i\phi}) \right)
```

# Arguments
- `j::Real`: The spin quantum number and can be a non-negative integer or half-integer
- `θ::Real`: rotation angle from z-axis
- `ϕ::Real`: rotation angle from x-axis

See also [`jmat`](@ref) and [`spin_state`](@ref).

# Reference
- [Robert Jones, Spin Coherent States and Statistical Physics](https://web.mit.edu/8.334/www/grades/projects/projects19/JonesRobert.pdf)
"""
function spin_coherent(j::Real, θ::Real, ϕ::Real)
    Sm = jmat(j, Val(:-))
    return exp(0.5 * θ * (Sm * exp(1im * ϕ) - Sm' * exp(-1im * ϕ))) * spin_state(j, j)
end

@doc raw"""
    bell_state(x::Int, z::Int)

Return the [Bell state](https://en.wikipedia.org/wiki/Bell_state) depending on the arguments `(x, z)`:
- `(0, 0)`: ``| \Phi^+ \rangle = ( |00\rangle + |11\rangle ) / \sqrt{2}``
- `(0, 1)`: ``| \Phi^- \rangle = ( |00\rangle - |11\rangle ) / \sqrt{2}``
- `(1, 0)`: ``| \Psi^+ \rangle = ( |01\rangle + |10\rangle ) / \sqrt{2}``
- `(1, 1)`: ``| \Psi^- \rangle = ( |01\rangle - |10\rangle ) / \sqrt{2}``

Here, `x = 1` (`z = 1`) means applying Pauli-``X`` ( Pauli-``Z``) unitary transformation on ``| \Phi^+ \rangle``.

# Example

```
julia> bell_state(0, 0)
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
bell_state(x::Int, z::Int) = bell_state(Val(x), Val(z))
bell_state(::Val{0}, ::Val{0}) = QuantumObject(ComplexF64[1, 0, 0, 1] / sqrt(2), Ket, [2, 2])
bell_state(::Val{0}, ::Val{1}) = QuantumObject(ComplexF64[1, 0, 0, -1] / sqrt(2), Ket, [2, 2])
bell_state(::Val{1}, ::Val{0}) = QuantumObject(ComplexF64[0, 1, 1, 0] / sqrt(2), Ket, [2, 2])
bell_state(::Val{1}, ::Val{1}) = QuantumObject(ComplexF64[0, 1, -1, 0] / sqrt(2), Ket, [2, 2])
bell_state(::Val{T1}, ::Val{T2}) where {T1,T2} = throw(ArgumentError("Invalid Bell state: $(T1), $(T2)"))

@doc raw"""
    singlet_state()

Return the two particle singlet state: ``\frac{1}{\sqrt{2}} ( |01\rangle - |10\rangle )``
"""
singlet_state() = QuantumObject(ComplexF64[0, 1, -1, 0] / sqrt(2), Ket, [2, 2])

@doc raw"""
    triplet_states()

Return a list of the two particle triplet states: 

- ``|11\rangle``
- ``( |01\rangle + |10\rangle ) / \sqrt{2}``
- ``|00\rangle``
"""
function triplet_states()
    return QuantumObject[
        QuantumObject(ComplexF64[0, 0, 0, 1], Ket, [2, 2]),
        QuantumObject(ComplexF64[0, 1, 1, 0] / sqrt(2), Ket, [2, 2]),
        QuantumObject(ComplexF64[1, 0, 0, 0], Ket, [2, 2]),
    ]
end

@doc raw"""
    w_state(n::Int)

Returns the `n`-qubit [W-state](https://en.wikipedia.org/wiki/W_state):

```math
\frac{1}{\sqrt{n}} \left( |100...0\rangle + |010...0\rangle + \cdots + |00...01\rangle \right)
```
"""
function w_state(n::Int)
    nzind = 2 .^ (0:(n-1)) .+ 1
    nzval = fill(ComplexF64(1 / sqrt(n)), n)
    return QuantumObject(SparseVector(2^n, nzind, nzval), Ket, fill(2, n))
end

@doc raw"""
    ghz_state(n::Int; d::Int=2)

Returns the generalized `n`-qudit [Greenberger–Horne–Zeilinger (GHZ) state](https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state):

```math
\frac{1}{\sqrt{d}} \sum_{i=0}^{d-1} | i \rangle \otimes \cdots \otimes | i \rangle
```

Here, `d` specifies the dimension of each qudit. Default to `d=2` (qubit).
"""
function ghz_state(n::Int; d::Int = 2)
    nzind = collect((0:(d-1)) .* Int((d^n - 1) / (d - 1)) .+ 1)
    nzval = ones(ComplexF64, d) / sqrt(d)
    return QuantumObject(SparseVector(d^n, nzind, nzval), Ket, fill(d, n))
end
