export dfd_mesolve, dsf_mesolve, dsf_mcsolve

### DYNAMICAL FOCK DIMENSION ###

function _reduce_dims(
    QO::AbstractArray{T},
    dims::Vector{<:Integer},
    sel::AbstractVector,
    reduce::AbstractVector,
) where {T}
    rd = dims
    nd = length(rd)
    rd_new = zero(rd)
    rd_new[sel] .= reduce
    @. rd_new = rd - rd_new

    if nd == 1
        ρmat = similar(QO, rd_new[1], rd_new[1])
        copyto!(ρmat, view(QO, 1:rd_new[1], 1:rd_new[1]))
    else
        ρmat = reshape(QO, reverse!(repeat(rd, 2))...)
        ρmat2 = similar(QO, reverse!(repeat(rd_new, 2))...)
        copyto!(ρmat2, view(ρmat, reverse!(repeat([1:n for n in rd_new], 2))...))
        ρmat = reshape(ρmat2, prod(rd_new), prod(rd_new))
    end

    return ρmat
end

function _increase_dims(
    QO::AbstractArray{T},
    dims::Vector{<:Integer},
    sel::AbstractVector,
    increase::AbstractVector,
) where {T}
    rd = dims
    nd = length(rd)
    rd_new = zero(rd)
    rd_new[sel] .= increase
    @. rd_new = rd + rd_new

    if nd == 1
        ρmat = similar(QO, rd_new[1], rd_new[1])
        fill!(selectdim(ρmat, 1, rd[1]+1:rd_new[1]), 0)
        fill!(selectdim(ρmat, 2, rd[1]+1:rd_new[1]), 0)
        copyto!(view(ρmat, 1:rd[1], 1:rd[1]), QO)
    else
        ρmat2 = similar(QO, reverse!(repeat(rd_new, 2))...)
        ρmat = reshape(QO, reverse!(repeat(rd, 2))...)
        for i in eachindex(sel)
            fill!(selectdim(ρmat2, nd - sel[i] + 1, rd[sel[i]]+1:rd_new[sel[i]]), 0)
            fill!(selectdim(ρmat2, 2 * nd - sel[i] + 1, rd[sel[i]]+1:rd_new[sel[i]]), 0)
        end
        copyto!(view(ρmat2, reverse!(repeat([1:n for n in rd], 2))...), ρmat)
        ρmat = reshape(ρmat2, prod(rd_new), prod(rd_new))
    end

    return ρmat
end

_dfd_set_pillow(dim)::Int = min(max(round(Int, 0.02 * dim), 1), 20)

function _DFDIncreaseReduceCondition(u, t, integrator)
    internal_params = integrator.p
    dim_list = internal_params.dim_list
    maxdims = internal_params.maxdims
    tol_list = internal_params.tol_list
    increase_list = internal_params.increase_list
    reduce_list = internal_params.reduce_list
    pillow_list = internal_params.pillow_list
    dfd_ρt_cache = internal_params.dfd_ρt_cache

    # I need this cache because I can't reshape directly the integrator.u
    copyto!(dfd_ρt_cache, u)

    @inbounds for i in eachindex(dim_list)
        maxdim_i = maxdims[i]
        dim_i = dim_list[i]
        pillow_i = pillow_list[i]
        if dim_i < maxdim_i && dim_i > 2 && maxdim_i != 0
            ρi = _ptrace_oper(vec2mat(dfd_ρt_cache), dim_list, [i])[1]
            @views res = norm(ρi[diagind(ρi)[end-pillow_i:end]], 1) * sqrt(dim_i) / pillow_i
            if res > tol_list[i]
                increase_list[i] = true
            elseif res < tol_list[i] * 1e-2 && dim_i > 3
                reduce_list[i] = true
            end
        end
    end
    return any(increase_list) || any(reduce_list)
end

function _DFDIncreaseReduceAffect!(integrator)
    internal_params = integrator.p
    H = internal_params.H_fun
    c_ops = internal_params.c_ops_fun
    e_ops = internal_params.e_ops_fun
    dim_list = internal_params.dim_list
    increase_list = internal_params.increase_list
    reduce_list = internal_params.reduce_list
    pillow_list = internal_params.pillow_list
    dim_list_evo_times = internal_params.dim_list_evo_times
    dim_list_evo = internal_params.dim_list_evo
    dfd_ρt_cache = internal_params.dfd_ρt_cache
    dfd_params = internal_params.dfd_params

    ρt = vec2mat(dfd_ρt_cache)

    @views pillow_increase = pillow_list[increase_list]
    @views pillow_reduce = pillow_list[reduce_list]
    dim_increase = findall(increase_list)
    dim_reduce = findall(reduce_list)

    if length(dim_increase) > 0
        ρt = _increase_dims(ρt, dim_list, dim_increase, pillow_increase)
        dim_list[dim_increase] .+= pillow_increase
    end
    if length(dim_reduce) > 0
        ρt = _reduce_dims(ρt, dim_list, dim_reduce, pillow_reduce)
        dim_list[dim_reduce] .-= pillow_reduce
    end

    @. pillow_list = _dfd_set_pillow(dim_list)
    fill!(increase_list, false)
    fill!(reduce_list, false)
    push!(dim_list_evo_times, integrator.t)
    push!(dim_list_evo, dim_list)

    e_ops2 = map(op -> mat2vec(get_data(op)'), e_ops(dim_list, dfd_params))
    L = liouvillian(H(dim_list, dfd_params), c_ops(dim_list, dfd_params)).data

    resize!(integrator, size(L, 1))
    copyto!(integrator.u, mat2vec(ρt))
    return integrator.p = merge(internal_params, (L = L, e_ops = e_ops2, dfd_ρt_cache = similar(integrator.u)))
end

function dfd_mesolveProblem(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T1},StateOpType},
    t_l::AbstractVector,
    c_ops::Function,
    maxdims::Vector{T2},
    dfd_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (dim_list) -> Vector{Vector{T1}}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    tol_list::Vector{<:Number} = fill(1e-8, length(maxdims)),
    kwargs...,
) where {T1,T2<:Integer,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    length(ψ0.dims) != length(maxdims) &&
        throw(DimensionMismatch("'dim_list' and 'maxdims' do not have the same dimension."))

    dim_list = deepcopy(ψ0.dims)
    H₀ = H(dim_list, dfd_params)
    c_ops₀ = c_ops(dim_list, dfd_params)
    e_ops₀ = e_ops(dim_list, dfd_params)

    dim_list_evo_times = [0.0]
    dim_list_evo = [dim_list]
    reduce_list = zeros(Bool, length(dim_list))
    increase_list = zeros(Bool, length(dim_list))
    pillow_list = _dfd_set_pillow.(dim_list)

    params2 = merge(
        params,
        (
            H_fun = H,
            c_ops_fun = c_ops,
            e_ops_fun = e_ops,
            dim_list = dim_list,
            maxdims = maxdims,
            tol_list = tol_list,
            reduce_list = reduce_list,
            increase_list = increase_list,
            pillow_list = pillow_list,
            dim_list_evo_times = dim_list_evo_times,
            dim_list_evo = dim_list_evo,
            dfd_ρt_cache = similar(ψ0.data, length(ψ0.data)^2),
            dfd_params = dfd_params,
        ),
    )

    cb_dfd = DiscreteCallback(_DFDIncreaseReduceCondition, _DFDIncreaseReduceAffect!, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb_dfd, kwargs2.callback),)) :
        merge(kwargs2, (callback = cb_dfd,))

    return mesolveProblem(H₀, ψ0, t_l, c_ops₀; e_ops = e_ops₀, alg = alg, H_t = H_t, params = params2, kwargs2...)
end

"""
    function dfd_mesolve(H::Function, ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function, maxdims::AbstractVector,
        dfd_params::NamedTuple=NamedTuple();
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Function=(dim_list) -> Vector{Vector{T1}}([]),
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        tol_list::Vector{<:Number}=fill(1e-8, length(maxdims)),
        kwargs...)

Time evolution of an open quantum system using master equation, dynamically changing the dimension of the Hilbert subspaces.
"""
function dfd_mesolve(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T1},StateOpType},
    t_l::AbstractVector,
    c_ops::Function,
    maxdims::Vector{T2},
    dfd_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (dim_list) -> Vector{Vector{T1}}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    tol_list::Vector{<:Number} = fill(1e-8, length(maxdims)),
    kwargs...,
) where {T1,T2<:Integer,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    dfd_prob = dfd_mesolveProblem(
        H,
        ψ0,
        t_l,
        c_ops,
        maxdims,
        dfd_params;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        tol_list = tol_list,
        kwargs...,
    )

    sol = solve(dfd_prob, alg)

    ρt = map(
        i -> QuantumObject(
            vec2mat(sol.u[i]),
            dims = sol.prob.p.dim_list_evo[searchsortedlast(sol.prob.p.dim_list_evo_times, sol.t[i])],
        ),
        eachindex(sol.t),
    )

    return TimeEvolutionSol(sol.t, ρt, sol.prob.p.expvals)
end

# Dynamical Shifted Fock mesolve

function _DSF_mesolve_Condition(u, t, integrator)
    internal_params = integrator.p
    op_l_vec = internal_params.op_l_vec
    δα_list = internal_params.δα_list

    condition = false
    @inbounds for i in eachindex(δα_list)
        op_vec = op_l_vec[i]
        δα = δα_list[i]
        Δα = dot(op_vec, u)
        if δα < abs(Δα)
            condition = true
        end
    end
    return condition
end

function _DSF_mesolve_Affect!(integrator)
    internal_params = integrator.p
    op_l = internal_params.op_l
    op_l_vec = internal_params.op_l_vec
    αt_list = internal_params.αt_list
    δα_list = internal_params.δα_list
    H = internal_params.H_fun
    c_ops = internal_params.c_ops_fun
    e_ops = internal_params.e_ops_fun
    e_ops_vec = internal_params.e_ops
    dsf_cache = internal_params.dsf_cache
    dsf_params = internal_params.dsf_params
    expv_cache = internal_params.expv_cache
    dsf_identity = internal_params.dsf_identity
    dsf_displace_cache_full = internal_params.dsf_displace_cache_full

    op_l_length = length(op_l)
    fill!(dsf_displace_cache_full.coefficients, 0)

    for i in eachindex(op_l)
        # op = op_l[i]
        op_vec = op_l_vec[i]
        αt = αt_list[i]
        δα = δα_list[i]
        Δα = dot(op_vec, integrator.u)

        if δα < abs(Δα)
            # Dᵢ = exp(Δα*op' - conj(Δα)*op)
            # copyto!(dsf_cache, integrator.u)
            # mul!(integrator.u, sprepost(Dᵢ', Dᵢ).data, dsf_cache)

            # This is equivalent to the code above, assuming that transpose(op) = adjoint(op)
            # Aᵢ = kron(Δα * op.data - conj(Δα) * op.data', dsf_identity) + kron(dsf_identity, conj(Δα) * op.data - Δα * op.data')

            # @. dsf_displace_cache_full[i] = Δα * dsf_displace_cache_left[i] - conj(Δα) * dsf_displace_cache_left_dag[i] + conj(Δα) * dsf_displace_cache_right[i] - Δα * dsf_displace_cache_right_dag[i]
            # Aᵢ = dsf_displace_cache_full[i]

            # dsf_cache .= integrator.u
            # arnoldi!(expv_cache, Aᵢ, dsf_cache)
            # expv!(integrator.u, expv_cache, one(αt), dsf_cache)

            dsf_displace_cache_full.coefficients[i] = Δα
            dsf_displace_cache_full.coefficients[i+op_l_length] = -conj(Δα)
            dsf_displace_cache_full.coefficients[i+2*op_l_length] = conj(Δα)
            dsf_displace_cache_full.coefficients[i+3*op_l_length] = -Δα

            αt_list[i] += Δα
        end
    end

    copyto!(dsf_cache, integrator.u)
    arnoldi!(expv_cache, dsf_displace_cache_full, dsf_cache)
    expv!(integrator.u, expv_cache, 1, dsf_cache)

    op_l2 = op_l .+ αt_list
    e_ops2 = e_ops(op_l2, dsf_params)
    _mat2vec_data = op -> mat2vec(get_data(op)')
    @. e_ops_vec = _mat2vec_data(e_ops2)
    return copyto!(internal_params.L, liouvillian(H(op_l2, dsf_params), c_ops(op_l2, dsf_params), dsf_identity).data)
end

function dsf_mesolveProblem(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T},StateOpType},
    t_l::AbstractVector,
    c_ops::Function,
    op_list::Vector{TOl},
    α0_l::Vector{<:Number} = zeros(length(op_list)),
    dsf_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (op_list, p) -> Vector{TOl}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
    krylov_dim::Int = min(5, cld(length(ket2dm(ψ0).data), 3)),
    kwargs...,
) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},TOl}
    op_l = op_list
    H₀ = H(op_l .+ α0_l, dsf_params)
    c_ops₀ = c_ops(op_l .+ α0_l, dsf_params)
    e_ops₀ = e_ops(op_l .+ α0_l, dsf_params)

    αt_list = convert(Vector{T}, α0_l)
    op_l_vec = map(op -> mat2vec(get_data(op)'), op_l)
    # Create the Krylov subspace with kron(H₀.data, H₀.data) just for initialize
    expv_cache = arnoldi(kron(H₀.data, H₀.data), mat2vec(ket2dm(ψ0).data), krylov_dim)
    dsf_identity = I(prod(H₀.dims))
    dsf_displace_cache_left = map(op -> Qobj(kron(op.data, dsf_identity)), op_l)
    dsf_displace_cache_left_dag = map(op -> Qobj(kron(sparse(op.data'), dsf_identity)), op_l)
    dsf_displace_cache_right = map(op -> Qobj(kron(dsf_identity, op.data)), op_l)
    dsf_displace_cache_right_dag = map(op -> Qobj(kron(dsf_identity, sparse(op.data'))), op_l)
    dsf_displace_cache_full = OperatorSum(
        zeros(length(op_l) * 4),
        vcat(
            dsf_displace_cache_left,
            dsf_displace_cache_left_dag,
            dsf_displace_cache_right,
            dsf_displace_cache_right_dag,
        ),
    )

    params2 = params
    params2 = merge(
        params,
        (
            H_fun = H,
            c_ops_fun = c_ops,
            e_ops_fun = e_ops,
            op_l = op_l,
            op_l_vec = op_l_vec,
            αt_list = αt_list,
            δα_list = δα_list,
            dsf_cache = similar(ψ0.data, length(ψ0.data)^2),
            expv_cache = expv_cache,
            dsf_identity = dsf_identity,
            dsf_params = dsf_params,
            dsf_displace_cache_full = dsf_displace_cache_full,
        ),
    )

    cb_dsf = DiscreteCallback(_DSF_mesolve_Condition, _DSF_mesolve_Affect!, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb_dsf, kwargs2.callback),)) :
        merge(kwargs2, (callback = cb_dsf,))

    return mesolveProblem(H₀, ψ0, t_l, c_ops₀; e_ops = e_ops₀, alg = alg, H_t = H_t, params = params2, kwargs2...)
end

"""
    function dsf_mesolve(H::Function,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function,
        op_list::Vector{TOl},
        α0_l::Vector{<:Number}=zeros(length(op_list)),
        dsf_params::NamedTuple=NamedTuple();
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Function=(op_list,p) -> Vector{TOl}([]),
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        δα_list::Vector{<:Number}=fill(0.2, length(op_list)),
        krylov_dim::Int=min(5,cld(length(ket2dm(ψ0).data), 3)),
        kwargs...)

Time evolution of an open quantum system using master equation and the Dynamical Shifted Fock algorithm.
"""
function dsf_mesolve(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T},StateOpType},
    t_l::AbstractVector,
    c_ops::Function,
    op_list::Vector{TOl},
    α0_l::Vector{<:Number} = zeros(length(op_list)),
    dsf_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (op_list, p) -> Vector{TOl}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
    krylov_dim::Int = min(5, cld(length(ket2dm(ψ0).data), 3)),
    kwargs...,
) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},TOl}
    dsf_prob = dsf_mesolveProblem(
        H,
        ψ0,
        t_l,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        δα_list = δα_list,
        krylov_dim = krylov_dim,
        kwargs...,
    )

    return mesolve(dsf_prob; alg = alg, kwargs...)
end

function dsf_mesolve(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T},StateOpType},
    t_l::AbstractVector,
    op_list::Vector{TOl},
    α0_l::Vector{<:Number} = zeros(length(op_list)),
    dsf_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (op_list, p) -> Vector{TOl}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
    krylov_dim::Int = min(5, cld(length(ket2dm(ψ0).data), 3)),
    kwargs...,
) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},TOl}
    c_ops = op_list -> Vector{TOl}([])
    return dsf_mesolve(
        H,
        ψ0,
        t_l,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        δα_list = δα_list,
        krylov_dim = krylov_dim,
        kwargs...,
    )
end

# Dynamical Shifted Fock mcsolve

function _DSF_mcsolve_Condition(u, t, integrator)
    internal_params = integrator.p
    op_l = internal_params.op_l
    δα_list = internal_params.δα_list

    ψt = u

    condition = false
    @inbounds for i in eachindex(op_l)
        op = op_l[i]
        δα = δα_list[i]
        Δα = dot(ψt, op.data, ψt) / dot(ψt, ψt)
        if δα < abs(Δα)
            condition = true
        end
    end
    return condition
end

function _DSF_mcsolve_Affect!(integrator)
    internal_params = integrator.p
    op_l = internal_params.op_l
    αt_list = internal_params.αt_list
    δα_list = internal_params.δα_list
    H = internal_params.H_fun
    c_ops = internal_params.c_ops_fun
    e_ops = internal_params.e_ops_fun
    e_ops0 = internal_params.e_ops_mc
    c_ops0 = internal_params.c_ops
    ψt = internal_params.dsf_cache1
    dsf_cache = internal_params.dsf_cache2
    expv_cache = internal_params.expv_cache
    dsf_params = internal_params.dsf_params
    dsf_displace_cache_full = internal_params.dsf_displace_cache_full

    copyto!(ψt, integrator.u)
    normalize!(ψt)

    op_l_length = length(op_l)
    fill!(dsf_displace_cache_full.coefficients, 0)

    for i in eachindex(op_l)
        op = op_l[i]
        αt = αt_list[i]
        δα = δα_list[i]
        Δα = dot(ψt, op.data, ψt)

        if δα < abs(Δα)
            # Dᵢ = exp(Δα*op' - conj(Δα)*op)
            # dsf_cache .= integrator.u
            # mul!(integrator.u, Dᵢ.data', dsf_cache)

            # Aᵢ = -Δα*op.data' + conj(Δα)*op.data
            # dsf_cache .= integrator.u
            # arnoldi!(expv_cache, Aᵢ, dsf_cache)
            # expv!(integrator.u, expv_cache, one(αt), dsf_cache)

            dsf_displace_cache_full.coefficients[i] = conj(Δα)
            dsf_displace_cache_full.coefficients[i+op_l_length] = -Δα

            αt_list[i] += Δα
        end
    end

    copyto!(dsf_cache, integrator.u)
    arnoldi!(expv_cache, dsf_displace_cache_full, dsf_cache)
    expv!(integrator.u, expv_cache, 1, dsf_cache)

    op_l2 = op_l .+ αt_list
    e_ops2 = e_ops(op_l2, dsf_params)
    c_ops2 = c_ops(op_l2, dsf_params)
    @. e_ops0 = get_data(e_ops2)
    @. c_ops0 = get_data(c_ops2)
    H_eff = H(op_l2, dsf_params).data - lmul!(convert(eltype(ψt), 0.5im), mapreduce(op -> op' * op, +, c_ops0))
    return mul!(internal_params.U, -1im, H_eff)
end

function _dsf_mcsolve_prob_func(prob, i, repeat)
    internal_params = prob.p

    prm = merge(
        internal_params,
        (
            U = copy(internal_params.U),
            e_ops_mc = copy(internal_params.e_ops_mc),
            c_ops = copy(internal_params.c_ops),
            expvals = similar(internal_params.expvals),
            cache_mc = similar(internal_params.cache_mc),
            weights_mc = similar(internal_params.weights_mc),
            cumsum_weights_mc = similar(internal_params.weights_mc),
            random_n = Ref(rand()),
            progr_mc = ProgressBar(size(internal_params.expvals, 2), enable = false),
            jump_times_which_idx = Ref(1),
            jump_times = similar(internal_params.jump_times),
            jump_which = similar(internal_params.jump_which),
            αt_list = copy(internal_params.αt_list),
            dsf_cache1 = similar(internal_params.dsf_cache1),
            dsf_cache2 = similar(internal_params.dsf_cache2),
            expv_cache = copy(internal_params.expv_cache),
            dsf_displace_cache_full = OperatorSum(
                copy(internal_params.dsf_displace_cache_full.coefficients),
                internal_params.dsf_displace_cache_full.operators,
            ),
        ),
    )

    return remake(prob, p = prm)
end

function dsf_mcsolveEnsembleProblem(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
    t_l::AbstractVector,
    c_ops::Function,
    op_list::Vector{TOl},
    α0_l::Vector{<:Number} = zeros(length(op_list)),
    dsf_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (op_list, p) -> Vector{TOl}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    krylov_dim::Int = min(5, cld(length(ψ0.data), 3)),
    kwargs...,
) where {T,TOl,TJC<:LindbladJumpCallbackType}
    op_l = op_list
    H₀ = H(op_l .+ α0_l, dsf_params)
    c_ops₀ = c_ops(op_l .+ α0_l, dsf_params)
    e_ops₀ = e_ops(op_l .+ α0_l, dsf_params)

    αt_list = convert(Vector{T}, α0_l)
    expv_cache = arnoldi(H₀.data, ψ0.data, krylov_dim)

    dsf_displace_cache = map(op -> Qobj(op.data), op_l)
    dsf_displace_cache_dag = map(op -> Qobj(sparse(op.data')), op_l)
    dsf_displace_cache_full = OperatorSum(zeros(length(op_l) * 2), vcat(dsf_displace_cache, dsf_displace_cache_dag))

    params2 = merge(
        params,
        (
            H_fun = H,
            c_ops_fun = c_ops,
            e_ops_fun = e_ops,
            op_l = op_l,
            αt_list = αt_list,
            δα_list = δα_list,
            dsf_cache1 = similar(ψ0.data),
            dsf_cache2 = similar(ψ0.data),
            expv_cache = expv_cache,
            dsf_params = dsf_params,
            dsf_displace_cache_full = dsf_displace_cache_full,
        ),
    )

    cb_dsf = DiscreteCallback(_DSF_mcsolve_Condition, _DSF_mcsolve_Affect!, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb_dsf, kwargs2.callback),)) :
        merge(kwargs2, (callback = cb_dsf,))

    return mcsolveEnsembleProblem(
        H₀,
        ψ0,
        t_l,
        c_ops₀;
        e_ops = e_ops₀,
        alg = alg,
        H_t = H_t,
        params = params2,
        jump_callback = jump_callback,
        prob_func = _dsf_mcsolve_prob_func,
        kwargs2...,
    )
end

"""
    function dsf_mcsolve(H::Function,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function,
        op_list::Vector{TOl},
        α0_l::Vector{<:Number}=zeros(length(op_list)),
        dsf_params::NamedTuple=NamedTuple();
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Function=(op_list,p) -> Vector{TOl}([]),
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        δα_list::Vector{<:Real}=fill(0.2, length(op_list)),
        n_traj::Int=1,
        ensemble_method=EnsembleThreads(),
        jump_callback::LindbladJumpCallbackType=ContinuousLindbladJumpCallback(),
        krylov_dim::Int=min(5,cld(length(ψ0.data), 3)),
        kwargs...)

Time evolution of a quantum system using the Monte Carlo wave function method
and the Dynamical Shifted Fock algorithm.
"""
function dsf_mcsolve(
    H::Function,
    ψ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
    t_l::AbstractVector,
    c_ops::Function,
    op_list::Vector{TOl},
    α0_l::Vector{<:Number} = zeros(length(op_list)),
    dsf_params::NamedTuple = NamedTuple();
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Function = (op_list, p) -> Vector{TOl}([]),
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
    n_traj::Int = 1,
    ensemble_method = EnsembleThreads(),
    jump_callback::TJC = ContinuousLindbladJumpCallback(),
    krylov_dim::Int = min(5, cld(length(ψ0.data), 3)),
    kwargs...,
) where {T,TOl,TJC<:LindbladJumpCallbackType}
    ens_prob_mc = dsf_mcsolveEnsembleProblem(
        H,
        ψ0,
        t_l,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        δα_list = δα_list,
        jump_callback = jump_callback,
        krylov_dim = krylov_dim,
        kwargs...,
    )

    return mcsolve(ens_prob_mc; alg = alg, n_traj = n_traj, ensemble_method = ensemble_method, kwargs...)
end
