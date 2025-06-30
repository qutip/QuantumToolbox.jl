#=
This file defines the energy restricted space structure.
=#

export EnrSpace, enr_state_dictionaries
export enr_identity, enr_fock, enr_destroy, enr_thermal_dm

@doc raw"""
    struct EnrSpace{N} <: AbstractSpace
        size::Int
        dims::NTuple{N,Int}
        n_excitations::Int
        state2idx::Dict{SVector{N,Int},Int}
        idx2state::Dict{Int,SVector{N,Int}}
    end

A structure that describes an excitation-number restricted (ENR) state space, where `N` is the number of sub-systems.

# Fields

- `size`: Number of states in the excitation-number restricted state space
- `dims`: A list of the number of states in each sub-system
- `n_excitations`: Maximum number of excitations
- `state2idx`: A dictionary for looking up a state index from a state (`SVector`)
- `idx2state`: A dictionary for looking up state (`SVector`) from a state index

# Example

To construct an `EnrSpace`, we only need to specify the `dims` and `n_excitations`, namely

```jldoctest
julia> dims = (2, 2, 3);

julia> n_excitations = 3;

julia> EnrSpace(dims, n_excitations)
EnrSpace((2, 2, 3), 3)
```
"""
struct EnrSpace{N} <: AbstractSpace
    size::Int
    dims::SVector{N,Int}
    n_excitations::Int
    state2idx::Dict{SVector{N,Int},Int}
    idx2state::Dict{Int,SVector{N,Int}}

    function EnrSpace(dims::Union{AbstractVector{T},NTuple{N,T}}, excitations::Int) where {T<:Integer,N}
        # all arguments will be checked in `enr_state_dictionaries`
        size, state2idx, idx2state = enr_state_dictionaries(dims, excitations)

        L = length(dims)
        return new{L}(size, SVector{L}(dims), excitations, state2idx, idx2state)
    end
end

function Base.show(io::IO, s::EnrSpace)
    print(io, "EnrSpace($(s.dims), $(s.n_excitations))")
    return nothing
end

Base.:(==)(s_enr1::EnrSpace, s_enr2::EnrSpace) = (s_enr1.size == s_enr2.size) && (s_enr1.dims == s_enr2.dims)

dimensions_to_dims(s_enr::EnrSpace) = s_enr.dims

@doc raw"""
    enr_state_dictionaries(dims, excitations)

Return the number of states, and lookup-dictionaries for translating a state (`SVector`) to a state index, and vice versa, for a system with a given number of components and maximum number of excitations.

# Arguments
- `dims::Union{AbstractVector,Tuple}`: A list of the number of states in each sub-system
- `excitations::Int`: Maximum number of excitations

# Returns
- `nstates`: Number of states
- `state2idx`: A dictionary for looking up a state index from a state (`SVector`)
- `idx2state`: A dictionary for looking up state (`SVector`) from a state index
"""
function enr_state_dictionaries(dims::Union{AbstractVector{T},NTuple{N,T}}, excitations::Int) where {T<:Integer,N}
    # argument checks
    _non_static_array_warning("dims", dims)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))
    all(>=(1), dims) || throw(DomainError(dims, "All the elements of dims must be non-zero integers (≥ 1)"))
    (excitations > 0) || throw(DomainError(excitations, "The argument excitations must be a non-zero integer (≥ 1)"))

    nvec = zeros(Int, L) # Vector
    nexc = 0

    # in the following, all `nvec` (Vector) will first be converted (copied) to SVector and then push to `result` 
    result = SVector{L,Int}[nvec]
    while true
        idx = L
        nvec[end] += 1
        nexc += 1
        (nvec[idx] < dims[idx]) && push!(result, nvec)
        while (nexc == excitations) || (nvec[idx] == dims[idx])
            idx -= 1

            # if idx < 1, break while-loop and return
            if idx < 1
                enr_size = length(result)
                return (enr_size, Dict(zip(result, 1:enr_size)), Dict(zip(1:enr_size, result)))
            end

            nexc -= nvec[idx+1] - 1
            nvec[idx+1] = 0
            nvec[idx] += 1
            (nvec[idx] < dims[idx]) && push!(result, nvec)
        end
    end
end

function enr_identity(dims::Union{AbstractVector{T},NTuple{N,T}}, excitations::Int) where {T<:Integer,N}
    s_enr = EnrSpace(dims, excitations)
    return enr_identity(s_enr)
end
enr_identity(s_enr::EnrSpace) = QuantumObject(Diagonal(ones(ComplexF64, s_enr.size)), Operator(), Dimensions(s_enr))

function enr_fock(
    dims::Union{AbstractVector{T},NTuple{N,T}},
    excitations::Int,
    state::AbstractVector{T};
    sparse::Union{Bool,Val} = Val(false),
) where {T<:Integer,N}
    s_enr = EnrSpace(dims, excitations)
    return enr_fock(s_enr, state, sparse = sparse)
end
function enr_fock(s_enr::EnrSpace, state::AbstractVector{T}; sparse::Union{Bool,Val} = Val(false)) where {T<:Integer}
    if getVal(sparse)
        array = sparsevec([s_enr.state2idx[[state...]]], [1.0 + 0im], s_enr.size)
    else
        j = s_enr.state2idx[state]
        array = [i == j ? 1.0 + 0im : 0.0 + 0im for i in 1:(s_enr.size)]
    end

    return QuantumObject(array, Ket(), s_enr)
end

function enr_destroy(dims::Union{AbstractVector{T},NTuple{N,T}}, excitations::Int) where {T<:Integer,N}
    s_enr = EnrSpace(dims, excitations)
    return enr_destroy(s_enr)
end
function enr_destroy(s_enr::EnrSpace{N}) where {N}
    D = s_enr.size
    idx2state = s_enr.idx2state
    state2idx = s_enr.state2idx

    a_ops = ntuple(i -> QuantumObject(spzeros(ComplexF64, D, D), Operator(), s_enr), N)

    for (n1, state1) in idx2state
        for (idx, s) in pairs(state1)
            # if s > 0, the annihilation operator of mode idx has a non-zero
            # entry with one less excitation in mode idx in the final state
            if s > 0
                state2 = Vector(state1)
                state2[idx] -= 1
                n2 = state2idx[state2]
                a_ops[idx][n2, n1] = √s
            end
        end
    end

    return a_ops
end

function enr_thermal_dm(
    dims::Union{AbstractVector{T1},NTuple{N,T1}},
    excitations::Int,
    n::Union{T2,AbstractVector{T2}},
) where {T1<:Integer,T2<:Real,N}
    s_enr = EnrSpace(dims, excitations)
    return enr_thermal_dm(s_enr, n)
end
function enr_thermal_dm(s_enr::EnrSpace{N}, n::Union{T,AbstractVector{T}}) where {N,T<:Real}
    if n isa Real
        nvec = fill(n, N)
    else
        (length(n) == N) || throw(ArgumentError("The length of the vector `n` should be the same as `dims`."))
        nvec = n
    end

    D = s_enr.size
    idx2state = s_enr.idx2state

    diags = ComplexF64[prod((nvec ./ (1 .+ nvec)) .^ idx2state[idx]) for idx in 1:D]

    diags /= sum(diags)

    return QuantumObject(Diagonal(diags), Operator(), s_enr)
end
