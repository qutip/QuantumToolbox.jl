#=
This file defines the energy restricted space structure.
=#

export EnrSpace, enr_state_dictionaries
export enr_identity, enr_fock, enr_destroy, enr_thermal_dm

struct EnrSpace{N} <: AbstractSpace
    size::Int
    dims::NTuple{N,Int}
    n_excitations::Int
    state2idx::Dict{SVector{N,Int},Int}
    idx2state::Dict{Int,SVector{N,Int}}

    function EnrSpace(dims::Union{Tuple,AbstractVector}, excitations::Int)
        _non_static_array_warning("dims", dims)
        dim_len = length(dims)
        dims_T = NTuple{dim_len}(dims)

        size, state2idx, idx2state = enr_state_dictionaries(dims, excitations)

        return new{dim_len}(size, dims_T, excitations, state2idx, idx2state)
    end
end

function Base.show(io::IO, s::EnrSpace)
    print(io, "EnrSpace($(s.dims), $(s.n_excitations))")
    return nothing
end

Base.:(==)(s_enr1::EnrSpace, s_enr2::EnrSpace) = (all([s_enr1.size, s_enr1.dims] .== [s_enr2.size, s_enr2.dims]))

dimensions_to_dims(s_enr::EnrSpace) = s_enr.dims

function enr_state_dictionaries(dims::Union{Tuple,AbstractVector}, excitations::Int)
    len = length(dims)
    nvec = zeros(Int, len)
    result = SVector{len,Int}[nvec] # in the following, all nvec will first be converted (copied) to SVector and then push to result 
    nexc = 0

    while true
        idx = len
        nvec[end] += 1
        nexc += 1
        if nvec[idx] < dims[idx]
            push!(result, nvec)
        end
        while (nexc == excitations) || (nvec[idx] == dims[idx])
            #nvec[idx] = 0
            idx -= 1
            if idx < 1
                enr_size = length(result)
                return (enr_size, Dict(zip(result, 1:enr_size)), Dict(zip(1:enr_size, result)))
            end

            nexc -= nvec[idx+1] - 1
            nvec[idx+1] = 0
            nvec[idx] += 1
            if nvec[idx] < dims[idx]
                push!(result, nvec)
            end
        end
    end
end

function enr_identity(dims::Union{Tuple,AbstractVector}, excitations::Int)
    s_enr = EnrSpace(dims, excitations)
    return QuantumObject(Diagonal(ones(ComplexF64, s_enr.size)), Operator(), Dimensions(s_enr))
end

function enr_fock(
    dims::Union{Tuple,AbstractVector},
    excitations::Int,
    state::AbstractVector;
    sparse::Union{Bool,Val} = Val(false),
)
    s_enr = EnrSpace(dims, excitations)
    if getVal(sparse)
        array = sparsevec([s_enr.state2idx[[state...]]], [1.0 + 0im], s_enr.size)
    else
        j = s_enr.state2idx[state]
        array = [i == j ? 1.0 + 0im : 0.0 + 0im for i in 1:(s_enr.size)]

        # s = zeros(ComplexF64, s_enr.size)
        # s[s_enr.state2idx[state]] += 1
        # s
    end

    return QuantumObject(array, Ket(), s_enr)
end

function enr_destroy(dims::Union{Tuple,AbstractVector}, excitations::Int)
    s_enr = EnrSpace(dims, excitations)
    N = s_enr.size
    idx2state = s_enr.idx2state
    state2idx = s_enr.state2idx

    a_ops = ntuple(i -> QuantumObject(spzeros(ComplexF64, N, N), Operator(), s_enr), length(dims))

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

function enr_thermal_dm(dims::Union{Tuple,AbstractVector}, excitations::Int, n::Union{Int,AbstractVector})
    if n isa Number
        nvec = Vector{typeof(n)}(n, length(dims))
    else
        length(n) == length(dims) || throw(ArgumentError("The Vector `n` has different length to `dims`."))
        nvec = n
    end

    s_enr = EnrSpace(dims, excitations)
    N = s_enr.size
    idx2state = s_enr.idx2state

    diags = [prod((nvec ./ (1 .+ nvec)) .^ idx2state[idx]) for idx in 1:N]

    diags /= sum(diags)

    return QuantumObject(Diagonal(diags), Operator(), s_enr)
end
