#=
Functions for generating (common) quantum states.
=#

export zero_ket, fock, basis, coherent
export fock_dm, coherent_dm, thermal_dm, maximally_mixed_dm, rand_dm
export spin_state, spin_coherent
export bell_state, singlet_state, triplet_states

@doc raw"""
    zero_ket(dimensions)

Returns a zero [`Ket`](@ref) vector with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.
"""
zero_ket(dimensions::Int) = QuantumObject(zeros(dimensions), Ket, [dimensions])
zero_ket(dimensions::Vector{Int}) = QuantumObject(zeros(prod(dimensions)), Ket, dimensions)

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

Generates a coherent state ``\ket{\alpha}``, which is defined as an eigenvector of the bosonic annihilation operator ``\hat{a} \ket{\alpha} = \alpha \ket{\alpha}``.
"""
function coherent(N::Int, α::T) where {T<:Number}
    a = destroy(N)
    return exp(α * a' - α' * a) * fock(N, 0)
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

Density matrix representation of a coherent state.

Constructed via outer product of [`coherent`](@ref).
"""
function coherent_dm(N::Int, α::T) where {T<:Number}
    ψ = coherent(N, α)
    return ψ * ψ'
end

@doc raw"""
    thermal_dm(N::Int, n::Real)

Density matrix for a thermal state (generating thermal state probabilities) with the following arguments:
- `N::Int`: Number of basis states in the Hilbert space
- `n::Real`: Expectation value for number of particles in the thermal state.
"""
function thermal_dm(N::Int, n::Real)
    β = log(1.0 / n + 1.0)
    N_list = Array{Float64}(0:N-1)
    data = exp.(-β .* N_list)
    return QuantumObject(spdiagm(0 => data ./ sum(data)), Operator, [N])
end

@doc raw"""
    maximally_mixed_dm(dimensions)

Returns the maximally mixed density matrix with given argument `dimensions`.

The `dimensions` can be either the following types:
- `dimensions::Int`: Number of basis states in the Hilbert space.
- `dimensions::Vector{Int}`: list of dimensions representing the each number of basis in the subsystems.
"""
maximally_mixed_dm(dimensions::Int) = QuantumObject(ones(dimensions, dimensions) / dimensions, Operator, [dimensions])
function maximally_mixed_dm(dimensions::Vector{Int})
    N = prod(dimensions)
    return QuantumObject(ones(N, N) / N, Operator, dimensions)
end

@doc raw"""
    rand_dm(N::Integer; dims::Vector{Int}=[N])

Generates a random density matrix ``\hat{\rho}``, with the property to be positive semi-definite and ``\textrm{Tr} \left[ \hat{\rho} \right] = 1``.

It is also possible to specify the list of dimensions `dims` if different subsystems are present.
"""
function rand_dm(N::Integer; dims::Vector{Int} = [N])
    ρ = rand(ComplexF64, N, N)
    ρ *= ρ'
    ρ /= tr(ρ)
    return QuantumObject(ρ; type = Operator, dims = dims)
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
    bell_state(state::String="00")

Return the corresponding Bell state:
- `"00"`: ``( |00\rangle + |11\rangle ) / \sqrt{2}``
- `"01"`: ``( |00\rangle - |11\rangle ) / \sqrt{2}``
- `"10"`: ``( |01\rangle + |10\rangle ) / \sqrt{2}``
- `"11"`: ``( |01\rangle - |10\rangle ) / \sqrt{2}``

# Example

```
julia> bell_state("00")
Quantum Object:   type=Ket   dims=[2, 2]   size=(4,)
4-element Vector{ComplexF64}:
 0.7071067811865475 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
function bell_state(state::String = "00")
    if state == "00"
        data = ComplexF64[1, 0, 0, 1]

    elseif state == "01"
        data = ComplexF64[1, 0, 0, -1]

    elseif state == "10"
        data = ComplexF64[0, 1, 1, 0]

    elseif state == "11"
        data = ComplexF64[0, 1, -1, 0]
    else
        throw(ArgumentError("Invalid state: $(state)"))
    end

    return QuantumObject(data / sqrt(2), Ket, [2, 2])
end

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
