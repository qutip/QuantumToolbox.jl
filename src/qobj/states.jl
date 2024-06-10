#=
Functions for generating (common) quantum states.
=#

export zero_ket, fock, basis, coherent
export fock_dm, coherent_dm, thermal_dm, maximally_mixed_dm, rand_dm

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
