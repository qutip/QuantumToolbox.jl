struct EnrSpace{N} <: AbstractSpace
    size::Int
    dims::SVector{N,Int}
    n_excitations::Int
    state2idx::Dict{}
    idx2state::Dict{}

    function EnrSpace(dims::Union{Tuple, AbstractVector}, excitations::Int)
        _non_static_array_warning("dims", dims)
        dim_len = length(dims)
        dims_T = NTuple{dim_len}(dims)

        size, state2idx, idx2state = enr_state_dictionaries(dims, excitations)
        

        return new{dim_len}(size, dims_T, excitations, state2idx, idx2state)
    end
end

function enr_state_dictionaries(dims::Union{Tuple, AbstractVector}, excitations::Int)
    len = length(dims)
    nvec = MVector{len}(zeros(Int, len))
    result = [copy(nvec)]
    nexc = 0

    while true
        idx = len
        nvec[end] += 1
        nexc += 1
        if nvec[idx] < dims[idx]
            push!(result, copy(nvec))
        end
        while (nexc == excitations) || (nvec[idx] == dims[idx])
            #nvec[idx] = 0
            idx -= 1
            if idx < 1
                enr_size = length(result)
                return (
                    enr_size,
                    Dict(zip(result, 1:enr_size)),
                    Dict(zip(1:enr_size, result))
                )
            end

            nexc -= nvec[idx+1] - 1
            nvec[idx+1] = 0
            nvec[idx] += 1
            if nvec[idx] < dims[idx]
                push!(result, copy(nvec))
            end
        end
    end

end

function enr_identity(dims::Union{Tuple, AbstractVector}, excitations::Int)
    s_enr = EnrSpace(dims, excitations)
    return QuantumObject(Diagonal(ones(ComplexF64, s_enr.size)), Operator(), Dimensions(s_enr))
end

function enr_fock(
        dims::Union{Tuple, AbstractVector}, excitations::Int, state::AbstractVector; 
        sparse::Union{Bool,Val} = Val(false)
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

function enr_destroy(dims::Union{Tuple, AbstractVector}, excitations::Int)
    s_enr = EnrSpace(dims, excitations)
    N = s_enr.size
    idx2state = s_enr.idx2state
    state2idx = s_enr.state2idx

    a_ops = [zeros(ComplexF64, N, N) for _ in 1:length(dims)]

    for (n1, state1) in idx2state
        for (idx, s) in pairs(state1)
            s > 0 || continue

            state2 = copy(state1)
            state2[idx] -= 1
            n2 = state2idx[state2]
            a_ops[idx][n2, n1] += √s
        end
    end

    return [QuantumObject(array, Operator(), s_enr) for array in a_ops]
end

function enr_thermal_dm(
        dims::Union{Tuple, AbstractVector}, 
        excitations::Int,
        n::Union{Int, AbstractVector}
    )
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
