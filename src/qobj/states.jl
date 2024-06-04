#=
Functions for generating (common) quantum states.
=#

export fock, basis, coherent
export rand_dm

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
"""
basis(N::Int, pos::Int = 0; dims::Vector{Int} = [N]) = fock(N, pos, dims = dims)

@doc raw"""
    coherent(N::Real, α::Number)

Generates a coherent state ``\ket{\alpha}``, which is defined as an eigenvector of the
bosonic annihilation operator ``\hat{a} \ket{\alpha} = \alpha \ket{\alpha}``.
"""
function coherent(N::Real, α::T) where {T<:Number}
    a = destroy(N)
    return exp(α * a' - α' * a) * fock(N, 0)
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
