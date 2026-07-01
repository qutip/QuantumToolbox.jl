#=
This file defines the excitation number restricted space structure.
=#

export EnrSpace, enr_state_dictionaries
export enr_fock, enr_thermal_dm, enr_destroy, enr_identity

@doc raw"""
    struct EnrSpace{N} <: AbstractSpace
        size::Int
        dims::NTuple{N,Int}
        n_excitations::Int
        excitation_weights::NTuple{N,Int}
        state2idx::Dict{SVector{N,Int},Int}
        idx2state::Dict{Int,SVector{N,Int}}
    end

A structure that describes an excitation number restricted (ENR) state space, where `N` is the number of sub-systems.

By default, the space contains all the states whose total number of excitations ``N = \sum_j n_j`` does not exceed `n_excitations`. More generally, one can restrict the space according to a weighted number operator ``N = \sum_j c_j n_j`` by passing the (positive integer) weights ``c_j`` through the keyword argument `excitation_weights`. This is useful whenever the Hamiltonian does not commute with the total number operator but conserves a weighted one, as it happens, e.g., for parametric down-conversion (``a^\dagger b^\dagger c + a b c^\dagger``), where a pump excitation in mode `c` is converted into one excitation in each of the modes `a` and `b`, so that ``n_a + n_b + 2 n_c`` is conserved.

# Fields

- `size`: Number of states in the excitation number restricted state space
- `dims`: A list of the number of states in each sub-system
- `n_excitations`: Maximum number of (weighted) excitations
- `excitation_weights`: The weights ``c_j`` associated to the number of excitations of each sub-system
- `state2idx`: A dictionary for looking up a state index from a state (`SVector`)
- `idx2state`: A dictionary for looking up state (`SVector`) from a state index

# Functions

With this `EnrSpace`, one can use the following functions to construct states or operators in the excitation number restricted (ENR) space:

- [`enr_fock`](@ref)
- [`enr_thermal_dm`](@ref)
- [`enr_destroy`](@ref)
- [`enr_identity`](@ref)

# Example

To construct an `EnrSpace`, we only need to specify the `dims` and `n_excitations`, namely

```jldoctest
julia> dims = (2, 2, 3);

julia> n_excitations = 3;

julia> EnrSpace(dims, n_excitations)
EnrSpace([2, 2, 3], 3)
```

One can also restrict the space according to a weighted number operator by specifying the `excitation_weights`, for example to keep only the states with ``n_1 + n_2 + 2 n_3 \leq 2``:

```jldoctest
julia> EnrSpace((2, 2, 2), 2; excitation_weights = (1, 1, 2))
EnrSpace([2, 2, 2], 2; excitation_weights = [1, 1, 2])
```

!!! note "Building number-conserving operators"
    The operators returned by [`enr_destroy`](@ref) act within the (truncated) restricted space. A product of such operators reproduces the corresponding full-space operator only if every intermediate state stays inside the restricted space. For a term that conserves the (weighted) number of excitations, this is guaranteed by writing it with the annihilation operators acting first (i.e. to the right) and obtaining the other half as its Hermitian conjugate. For instance, the parametric down-conversion interaction should be built as `g * (a' * b' * c + (a' * b' * c)')` rather than reordering the factors as `g * (a' * b' * c + a * b * c')`, since the latter would create an intermediate state outside the restricted space.

!!! warning "Beware of type-stability!"
    It is highly recommended to use `EnrSpace(dims, n_excitations)` with `dims` (and `excitation_weights`) as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
struct EnrSpace{N} <: AbstractSpace
    size::Int
    dims::SVector{N, Int}
    n_excitations::Int
    excitation_weights::SVector{N, Int}
    state2idx::Dict{SVector{N, Int}, Int}
    idx2state::Dict{Int, SVector{N, Int}}

    function EnrSpace(
            dims::AbstractVecOrTuple{T},
            n_excitations::Int;
            excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
        ) where {T <: Integer}
        # all arguments will be checked in `enr_state_dictionaries`
        size, state2idx, idx2state = enr_state_dictionaries(dims, n_excitations; excitation_weights = excitation_weights)

        L = length(dims)
        return new{L}(size, SVector{L}(dims), n_excitations, SVector{L}(excitation_weights), state2idx, idx2state)
    end
end

function Base.show(io::IO, s::EnrSpace)
    if all(==(1), s.excitation_weights)
        print(io, "EnrSpace($(s.dims), $(s.n_excitations))")
    else
        print(io, "EnrSpace($(s.dims), $(s.n_excitations); excitation_weights = $(s.excitation_weights))")
    end
    return nothing
end

Base.length(::EnrSpace{N}) where {N} = N

Base.:(==)(s_enr1::EnrSpace, s_enr2::EnrSpace) =
    (s_enr1.size == s_enr2.size) && (s_enr1.dims == s_enr2.dims) && (s_enr1.excitation_weights == s_enr2.excitation_weights)

dimensions_to_dims(s_enr::EnrSpace) = s_enr.dims

get_size(s_enr::EnrSpace) = s_enr.size

@doc raw"""
    enr_state_dictionaries(dims, n_excitations; excitation_weights = ntuple(_ -> 1, length(dims)))

Return the number of states, and lookup-dictionaries for translating a state (`SVector`) to a state index, and vice versa, for a system with a given number of components and maximum number of excitations.

The returned states are the ones whose weighted number of excitations ``\sum_j c_j n_j`` (with weights ``c_j`` given by `excitation_weights`) does not exceed `n_excitations`. By default all the weights are equal to one, recovering the standard total-excitation restriction ``\sum_j n_j \leq`` `n_excitations`.

# Arguments
- `dims::Union{AbstractVector,Tuple}`: A list of the number of states in each sub-system
- `n_excitations::Int`: Maximum number of (weighted) excitations
- `excitation_weights::Union{AbstractVector,Tuple}`: The (positive integer) weights ``c_j`` associated to the number of excitations of each sub-system

# Returns
- `nstates`: Number of states
- `state2idx`: A dictionary for looking up a state index from a state (`SVector`)
- `idx2state`: A dictionary for looking up state (`SVector`) from a state index
"""
function enr_state_dictionaries(
        dims::AbstractVecOrTuple{T},
        n_excitations::Int;
        excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
    ) where {T <: Integer}
    # argument checks
    _non_static_array_warning("dims", dims)
    _non_static_array_warning("excitation_weights", excitation_weights)
    L = length(dims)
    (L > 0) || throw(DomainError(dims, "The argument dims must be of non-zero length"))
    all(>=(1), dims) || throw(DomainError(dims, "All the elements of dims must be non-zero integers (≥ 1)"))
    (n_excitations > 0) ||
        throw(DomainError(n_excitations, "The argument n_excitations must be a non-zero integer (≥ 1)"))
    (length(excitation_weights) == L) ||
        throw(DimensionMismatch("The argument excitation_weights must have the same length as dims."))
    all(>=(1), excitation_weights) ||
        throw(DomainError(excitation_weights, "All the elements of excitation_weights must be non-zero integers (≥ 1)"))

    # Recursively enumerate, in lexicographic order (first sub-system being the most
    # significant), all the states `nvec` in the number basis such that
    #   0 ≤ nvec[j] ≤ dims[j] - 1   and   ∑ⱼ excitation_weights[j] * nvec[j] ≤ n_excitations
    result = SVector{L, Int}[]
    nvec = zeros(Int, L) # Vector
    _enr_append_states!(result, nvec, dims, excitation_weights, n_excitations, 1)

    enr_size = length(result)
    return (enr_size, Dict(zip(result, 1:enr_size)), Dict(zip(1:enr_size, result)))
end

# Append to `result` all the (weighted) excitation-restricted states obtained by choosing the
# occupation of sub-systems `idx, idx+1, …` given that `remaining` excitations are still available
# (see [`enr_state_dictionaries`](@ref)).
function _enr_append_states!(result, nvec::Vector{Int}, dims, excitation_weights, remaining::Int, idx::Int)
    if idx > length(nvec)
        push!(result, nvec) # converted (copied) to SVector when pushed to the typed `result`
        return nothing
    end

    # the occupation of sub-system `idx` is bounded both by its dimension and by the
    # remaining excitation budget divided by its weight
    n_max = min(dims[idx] - 1, remaining ÷ excitation_weights[idx])
    for n in 0:n_max
        nvec[idx] = n
        _enr_append_states!(result, nvec, dims, excitation_weights, remaining - excitation_weights[idx] * n, idx + 1)
    end
    return nothing
end

@doc raw"""
    enr_fock([T::Type=ComplexF64,] dims::Union{AbstractVector,Tuple}, n_excitations::Int, state::AbstractVector; sparse::Union{Bool,Val}=Val(false), excitation_weights::Union{AbstractVector,Tuple})
    enr_fock([T::Type=ComplexF64,] s_enr::EnrSpace, state::AbstractVector; sparse::Union{Bool,Val}=Val(false))

Generate the Fock state representation ([`Ket`](@ref)) in an excitation number restricted state space ([`EnrSpace`](@ref)) with element type `T = ComplexF64` (default).

The arguments `dims`, `n_excitations` and `excitation_weights` are used to generate [`EnrSpace`](@ref) (see its docstring for the meaning of `excitation_weights`).

The `state` argument is a list of integers that specifies the state (in the number basis representation) for which to generate the Fock state representation.

!!! warning "Beware of type-stability!"
    It is highly recommended to use `enr_fock(dims, n_excitations, state)` with `dims` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function enr_fock(
        ::Type{T},
        dims::AbstractVecOrTuple{Td},
        n_excitations::Int,
        state::AbstractVector{Td};
        sparse::Union{Bool, Val} = Val(false),
        excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
    ) where {T <: Number, Td <: Integer}
    s_enr = EnrSpace(dims, n_excitations; excitation_weights = excitation_weights)
    return enr_fock(T, s_enr, state; sparse)
end
function enr_fock(::Type{T}, s_enr::EnrSpace, state::AbstractVector{Td}; sparse::Union{Bool, Val} = Val(false)) where {T <: Number, Td <: Integer}
    if getVal(sparse)
        array = sparsevec([s_enr.state2idx[[state...]]], [one(T)], s_enr.size)
    else
        j = s_enr.state2idx[state]
        z0 = zero(T)
        array = [i == j ? one(T) : z0 for i in 1:(s_enr.size)]
    end

    return QuantumObject(array, Ket(), s_enr)
end
enr_fock(
    dims::AbstractVecOrTuple{Td},
    n_excitations::Int,
    state::AbstractVector{Td};
    sparse::Union{Bool, Val} = Val(false),
    excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
) where {Td <: Integer} =
    enr_fock(ComplexF64, dims, n_excitations, state; sparse, excitation_weights)
enr_fock(s_enr::EnrSpace, state::AbstractVector{Td}; sparse::Union{Bool, Val} = Val(false)) where {Td <: Integer} = enr_fock(ComplexF64, s_enr, state; sparse)

@doc raw"""
    enr_thermal_dm(dims::Union{AbstractVector,Tuple}, n_excitations::Int, n::Union{Real,AbstractVector}; sparse::Union{Bool,Val}=Val(false), excitation_weights::Union{AbstractVector,Tuple})
    enr_thermal_dm(s_enr::EnrSpace, n::Union{Real,AbstractVector}; sparse::Union{Bool,Val}=Val(false))

Generate the thermal state (density [`Operator`](@ref)) in an excitation number restricted state space ([`EnrSpace`](@ref)) with element type same as `n`.

The arguments `dims`, `n_excitations` and `excitation_weights` are used to generate [`EnrSpace`](@ref) (see its docstring for the meaning of `excitation_weights`).

The argument `n` is a list that specifies the expectation values for number of particles in each sub-system. If `n` is specified as a real number, it will apply to each sub-system.

!!! warning "Beware of type-stability!"
    It is highly recommended to use `enr_thermal_dm(dims, n_excitations, n)` with `dims` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function enr_thermal_dm(
        dims::AbstractVecOrTuple{T1},
        n_excitations::Int,
        n::Union{T2, AbstractVector{T2}};
        sparse::Union{Bool, Val} = Val(false),
        excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
    ) where {T1 <: Integer, T2 <: Real}
    s_enr = EnrSpace(dims, n_excitations; excitation_weights = excitation_weights)
    return enr_thermal_dm(s_enr, n; sparse)
end
function enr_thermal_dm(
        s_enr::EnrSpace{N},
        n::Union{T, AbstractVector{T}};
        sparse::Union{Bool, Val} = Val(false),
    ) where {N, T <: Real}
    if n isa Real
        nvec = fill(n, N)
    else
        (length(n) == N) || throw(ArgumentError("The length of the vector `n` should be the same as `dims`."))
        nvec = n
    end

    D = s_enr.size
    idx2state = s_enr.idx2state

    β = @. log(1 + 1 / nvec) # here makes element type become float
    P = [
        prod(_Boltzmann_weight(β[k], n_excite) for (k, n_excite) in pairs(idx2state[idx])) for idx in 1:D
    ]
    P /= sum(P)
    if getVal(sparse)
        return QuantumObject(spdiagm(0 => P), Operator(), s_enr)
    else
        return QuantumObject(diagm(0 => P), Operator(), s_enr)
    end
end

@doc raw"""
    enr_destroy([T::Type=ComplexF64,] dims::Union{AbstractVector,Tuple}, n_excitations::Int; excitation_weights::Union{AbstractVector,Tuple})
    enr_destroy([T::Type=ComplexF64,] s_enr::EnrSpace)

Generate a `Tuple` of annihilation operators for each sub-system in an excitation number restricted state space ([`EnrSpace`](@ref)) with element type `T = ComplexF64` (default). Thus, the return `Tuple` will have the same length as `dims`.

The arguments `dims`, `n_excitations` and `excitation_weights` are used to generate [`EnrSpace`](@ref) (see its docstring for the meaning of `excitation_weights`).

!!! warning "Beware of type-stability!"
    It is highly recommended to use `enr_destroy(dims, n_excitations)` with `dims` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function enr_destroy(
        ::Type{T},
        dims::AbstractVecOrTuple{Td},
        n_excitations::Int;
        excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
    ) where {T <: FloatOrComplex, Td <: Integer}
    s_enr = EnrSpace(dims, n_excitations; excitation_weights = excitation_weights)
    return enr_destroy(T, s_enr)
end
function enr_destroy(::Type{T}, s_enr::EnrSpace{N}) where {T <: FloatOrComplex, N}
    D = s_enr.size
    idx2state = s_enr.idx2state
    state2idx = s_enr.state2idx

    I_list = [Int64[] for _ in 1:N]
    J_list = [Int64[] for _ in 1:N]
    V_list = [T[] for _ in 1:N]

    for (n1, state1) in idx2state
        for (idx, s) in pairs(state1)
            # if s > 0, the annihilation operator of mode idx has a non-zero
            # entry with one less excitation in mode idx in the final state
            if s > 0
                state2 = Vector(state1)
                state2[idx] -= 1
                n2 = state2idx[state2]
                push!(I_list[idx], n2)
                push!(J_list[idx], n1)
                push!(V_list[idx], sqrt(T(s)))
            end
        end
    end

    return ntuple(i -> QuantumObject(sparse(I_list[i], J_list[i], V_list[i], D, D), Operator(), s_enr), Val(N))
end
enr_destroy(
    dims::AbstractVecOrTuple{Td},
    n_excitations::Int;
    excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
) where {Td <: Integer} = enr_destroy(ComplexF64, dims, n_excitations; excitation_weights)
enr_destroy(s_enr::EnrSpace{N}) where {N} = enr_destroy(ComplexF64, s_enr)

@doc raw"""
    enr_identity([T::Type=ComplexF64,] dims::Union{AbstractVector,Tuple}, n_excitations::Int; excitation_weights::Union{AbstractVector,Tuple})
    enr_identity([T::Type=ComplexF64,] s_enr::EnrSpace)

Generate the identity operator in an excitation number restricted state space ([`EnrSpace`](@ref)) with element type `T = ComplexF64` (default).

The arguments `dims`, `n_excitations` and `excitation_weights` are used to generate [`EnrSpace`](@ref) (see its docstring for the meaning of `excitation_weights`).

!!! warning "Beware of type-stability!"
    It is highly recommended to use `enr_identity(dims, n_excitations)` with `dims` as `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to keep type stability. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function enr_identity(
        ::Type{T},
        dims::AbstractVecOrTuple{Td},
        n_excitations::Int;
        excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
    ) where {T <: Number, Td <: Integer}
    s_enr = EnrSpace(dims, n_excitations; excitation_weights = excitation_weights)
    return enr_identity(T, s_enr)
end
enr_identity(::Type{T}, s_enr::EnrSpace) where {T <: Number} = QuantumObject(Diagonal(ones(T, s_enr.size)), Operator(), s_enr)
enr_identity(
    dims::AbstractVecOrTuple{Td},
    n_excitations::Int;
    excitation_weights::AbstractVecOrTuple = ntuple(_ -> 1, length(dims)),
) where {Td <: Integer} = enr_identity(ComplexF64, dims, n_excitations; excitation_weights)
enr_identity(s_enr::EnrSpace) = enr_identity(ComplexF64, s_enr)
