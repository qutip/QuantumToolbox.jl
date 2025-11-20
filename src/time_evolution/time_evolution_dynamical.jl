export dfd_mesolve, dsf_mesolve, dsf_mcsolve

### DYNAMICAL FOCK DIMENSION ###

function _reduce_dims(
        QO::AbstractArray{T},
        dims::Union{SVector{N, DT}, MVector{N, DT}},
        sel,
        reduce,
    ) where {T, N, DT <: Integer}
    n_d = length(dims)
    dims_new = zero(dims)
    dims_new[sel] .= reduce
    @. dims_new = dims - dims_new

    if n_d == 1
        ρmat = similar(QO, dims_new[1], dims_new[1])
        copyto!(ρmat, view(QO, 1:dims_new[1], 1:dims_new[1]))
    else
        ρmat = reshape(QO, reverse(vcat(dims, dims))...)
        ρmat2 = similar(QO, reverse(vcat(dims_new, dims_new))...)
        copyto!(ρmat2, view(ρmat, reverse!(repeat([1:n for n in dims_new], 2))...))
        ρmat = reshape(ρmat2, prod(dims_new), prod(dims_new))
    end

    return ρmat
end

function _increase_dims(
        QO::AbstractArray{T},
        dims::Union{SVector{N, DT}, MVector{N, DT}},
        sel,
        increase,
    ) where {T, N, DT <: Integer}
    n_d = length(dims)
    dims_new = MVector(zero(dims)) # Mutable SVector
    dims_new[sel] .= increase
    @. dims_new = dims + dims_new

    if n_d == 1
        ρmat = similar(QO, dims_new[1], dims_new[1])
        fill!(selectdim(ρmat, 1, (dims[1] + 1):dims_new[1]), 0)
        fill!(selectdim(ρmat, 2, (dims[1] + 1):dims_new[1]), 0)
        copyto!(view(ρmat, 1:dims[1], 1:dims[1]), QO)
    else
        ρmat2 = similar(QO, reverse(vcat(dims_new, dims_new))...)
        ρmat = reshape(QO, reverse(vcat(dims, dims))...)
        for i in eachindex(sel)
            fill!(selectdim(ρmat2, n_d - sel[i] + 1, (dims[sel[i]] + 1):dims_new[sel[i]]), 0)
            fill!(selectdim(ρmat2, 2 * n_d - sel[i] + 1, (dims[sel[i]] + 1):dims_new[sel[i]]), 0)
        end
        copyto!(view(ρmat2, reverse!(repeat([1:n for n in dims], 2))...), ρmat)
        ρmat = reshape(ρmat2, prod(dims_new), prod(dims_new))
    end

    return ρmat
end

_dfd_set_pillow(dim)::Int = min(max(round(Int, 0.02 * dim), 1), 20)

function _DFDIncreaseReduceCondition(u, t, integrator)
    params = integrator.p
    dim_list = params.dim_list
    maxdims = params.maxdims
    tol_list = params.tol_list
    increase_list = params.increase_list
    reduce_list = params.reduce_list
    pillow_list = params.pillow_list
    dfd_ρt_cache = params.dfd_ρt_cache

    # I need this cache because I can't reshape directly the integrator.u
    copyto!(dfd_ρt_cache, u)

    @inbounds for i in eachindex(dim_list)
        maxdim_i = maxdims[i]
        dim_i = dim_list[i]
        pillow_i = pillow_list[i]
        if dim_i < maxdim_i && dim_i > 2 && maxdim_i != 0
            ρi = _ptrace_oper(vec2mat(dfd_ρt_cache), dim_list, SVector(i))[1]
            @views res = norm(ρi[diagind(ρi)[(end - pillow_i):end]], 1) * sqrt(dim_i) / pillow_i
            if res > tol_list[i]
                increase_list[i] = true
            elseif res < tol_list[i] * 1.0e-2 && dim_i > 3
                reduce_list[i] = true
            end
        end
    end
    return any(increase_list) || any(reduce_list)
end

function _DFDIncreaseReduceAffect!(integrator)
    params = integrator.p
    H = params.H_fun
    c_ops = params.c_ops_fun
    e_ops = params.e_ops_fun
    dim_list = params.dim_list
    increase_list = params.increase_list
    reduce_list = params.reduce_list
    pillow_list = params.pillow_list
    dim_list_evo_times = params.dim_list_evo_times
    dim_list_evo = params.dim_list_evo
    dfd_ρt_cache = params.dfd_ρt_cache
    dfd_params = params.dfd_params

    ρt = vec2mat(dfd_ρt_cache)

    pillow_increase = pillow_list[increase_list] # TODO: This returns a Vector. Find a way to return an SVector or NTuple
    pillow_reduce = pillow_list[reduce_list]
    dim_increase = findall(increase_list) # TODO: This returns a Vector. Find a way to return an SVector or NTuple
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
    push!(dim_list_evo, copy(dim_list))

    e_ops2 = map(op -> mat2vec(get_data(op)'), e_ops(dim_list, dfd_params))
    L = liouvillian(H(dim_list, dfd_params), c_ops(dim_list, dfd_params)).data

    resize!(integrator, size(L, 1))
    copyto!(integrator.u, mat2vec(ρt))
    # By doing this, we are assuming that the system is time-independent and f is a MatrixOperator
    integrator.f = ODEFunction{true, FullSpecialize}(MatrixOperator(L))
    integrator.p = merge(params, (dfd_ρt_cache = similar(integrator.u),))
    _mesolve_callbacks_new_e_ops!(integrator, e_ops2)

    return nothing
end

function dfd_mesolveProblem(
        H::Function,
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        c_ops::Function,
        maxdims::Vector{T2},
        dfd_params::NamedTuple = NamedTuple();
        e_ops::Function = (dim_list) -> Vector{Vector{eltype(ψ0)}}([]),
        params::NamedTuple = NamedTuple(),
        tol_list::Vector{<:Number} = fill(1.0e-8, length(maxdims)),
        kwargs...,
    ) where {T2 <: Integer, StateOpType <: Union{Ket, Operator}}
    length(ψ0.dimensions) != length(maxdims) &&
        throw(DimensionMismatch("`dim_list` and `maxdims` do not have the same dimension."))

    dim_list = MVector(ψ0.dims)
    H₀ = H(dim_list, dfd_params)
    c_ops₀ = c_ops(dim_list, dfd_params)
    e_ops₀ = e_ops(dim_list, dfd_params)

    dim_list_evo_times = [0.0]
    dim_list_evo = [copy(dim_list)]
    reduce_list = MVector(ntuple(i -> false, Val(length(dim_list))))
    increase_list = MVector(ntuple(i -> false, Val(length(dim_list))))
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

    return mesolveProblem(H₀, ψ0, tlist, c_ops₀; e_ops = e_ops₀, params = params2, kwargs2...)
end

@doc raw"""
    dfd_mesolve(H::Function, ψ0::QuantumObject,
        tlist::AbstractVector, c_ops::Function, maxdims::AbstractVector,
        dfd_params::NamedTuple=NamedTuple();
        alg::AbstractODEAlgorithm=DP5(),
        e_ops::Function=(dim_list) -> Vector{Vector{T1}}([]),
        params::NamedTuple=NamedTuple(),
        tol_list::Vector{<:Number}=fill(1e-8, length(maxdims)),
        kwargs...)

Time evolution of an open quantum system using master equation, dynamically changing the dimension of the Hilbert subspaces.

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function dfd_mesolve(
        H::Function,
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        c_ops::Function,
        maxdims::Vector{T2},
        dfd_params::NamedTuple = NamedTuple();
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Function = (dim_list) -> Vector{Vector{eltype(ψ0)}}([]),
        params::NamedTuple = NamedTuple(),
        tol_list::Vector{<:Number} = fill(1.0e-8, length(maxdims)),
        kwargs...,
    ) where {T2 <: Integer, StateOpType <: Union{Ket, Operator}}
    dfd_prob = dfd_mesolveProblem(
        H,
        ψ0,
        tlist,
        c_ops,
        maxdims,
        dfd_params;
        e_ops = e_ops,
        params = params,
        tol_list = tol_list,
        kwargs...,
    )

    sol = solve(dfd_prob.prob, alg)

    ρt = map(eachindex(sol.t)) do i
        idx = findfirst(>=(sol.t[i]), sol.prob.p.dim_list_evo_times)
        idx2 = isnothing(idx) ? length(sol.prob.p.dim_list_evo) : (idx == 1 ? 1 : idx - 1)

        return QuantumObject(vec2mat(sol.u[i]), type = Operator(), dims = sol.prob.p.dim_list_evo[idx2])
    end

    return TimeEvolutionSol(
        dfd_prob.times,
        sol.t,
        ρt,
        _get_expvals(sol, SaveFuncMESolve),
        sol.retcode,
        sol.alg,
        sol.prob.kwargs[:abstol],
        sol.prob.kwargs[:reltol],
    )
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
    op_list = internal_params.op_list
    op_l_vec = internal_params.op_l_vec
    αt_list = internal_params.αt_list
    δα_list = internal_params.δα_list
    H = internal_params.H_fun
    c_ops = internal_params.c_ops_fun
    e_ops = internal_params.e_ops_fun
    dsf_cache = internal_params.dsf_cache
    dsf_params = internal_params.dsf_params
    expv_cache = internal_params.expv_cache
    dsf_displace_cache_full = internal_params.dsf_displace_cache_full

    op_l_length = length(op_list)

    for i in eachindex(op_list)
        # op = op_list[i]
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

            dsf_displace_cache_full.ops[i].λ.val = Δα
            dsf_displace_cache_full.ops[i + op_l_length].λ.val = -conj(Δα)
            dsf_displace_cache_full.ops[i + 2 * op_l_length].λ.val = conj(Δα)
            dsf_displace_cache_full.ops[i + 3 * op_l_length].λ.val = -Δα

            αt_list[i] += Δα
        else
            dsf_displace_cache_full.ops[i].λ.val = 0
            dsf_displace_cache_full.ops[i + op_l_length].λ.val = 0
            dsf_displace_cache_full.ops[i + 2 * op_l_length].λ.val = 0
            dsf_displace_cache_full.ops[i + 3 * op_l_length].λ.val = 0
        end
    end

    copyto!(dsf_cache, integrator.u)
    arnoldi!(expv_cache, dsf_displace_cache_full, dsf_cache)
    expv!(integrator.u, expv_cache, 1, dsf_cache)

    op_l2 = op_list .+ αt_list
    e_ops2 = e_ops(op_l2, dsf_params)
    _mesolve_callbacks_new_e_ops!(integrator, [_generate_mesolve_e_op(op) for op in e_ops2])

    # By doing this, we are assuming that all the arguments of ODEFunction are the default ones
    integrator.f =
        ODEFunction{true, FullSpecialize}(_mesolve_make_L_QobjEvo(H(op_l2, dsf_params), c_ops(op_l2, dsf_params)).data)
    return u_modified!(integrator, true)
end

function dsf_mesolveProblem(
        H::Function,
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        c_ops::Function,
        op_list::Union{AbstractVector, Tuple},
        α0_l::Vector{<:Number} = zeros(length(op_list)),
        dsf_params::NamedTuple = NamedTuple();
        e_ops::Function = (op_list, p) -> (),
        params::NamedTuple = NamedTuple(),
        δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
        krylov_dim::Int = max(6, min(10, cld(length(ket2dm(ψ0).data), 4))),
        kwargs...,
    ) where {StateOpType <: Union{Ket, Operator}}
    op_list = deepcopy(op_list)
    H₀ = H(op_list .+ α0_l, dsf_params)
    c_ops₀ = c_ops(op_list .+ α0_l, dsf_params)
    e_ops₀ = e_ops(op_list .+ α0_l, dsf_params)

    T = eltype(ψ0)

    αt_list = convert(Vector{T}, α0_l)
    op_l_vec = map(op -> mat2vec(get_data(op)'), op_list)

    dsf_identity = Eye(hilbert_dimensions_to_size(H₀.dimensions)[1])

    # Create the Krylov subspace just for initialize
    expv_cache = arnoldi(kron(dsf_identity, dsf_identity), mat2vec(ket2dm(ψ0).data), krylov_dim)

    dsf_displace_cache_left = sum(op -> ScalarOperator(one(T)) * MatrixOperator(kron(op.data, dsf_identity)), op_list)
    dsf_displace_cache_left_dag =
        sum(op -> ScalarOperator(one(T)) * MatrixOperator(kron(sparse(op.data'), dsf_identity)), op_list)
    dsf_displace_cache_right = sum(op -> ScalarOperator(one(T)) * MatrixOperator(kron(dsf_identity, op.data)), op_list)
    dsf_displace_cache_right_dag =
        sum(op -> ScalarOperator(one(T)) * MatrixOperator(kron(dsf_identity, sparse(op.data'))), op_list)
    dsf_displace_cache_full =
        dsf_displace_cache_left + dsf_displace_cache_left_dag + dsf_displace_cache_right + dsf_displace_cache_right_dag

    params2 = merge(
        params,
        (
            H_fun = H,
            c_ops_fun = c_ops,
            e_ops_fun = e_ops,
            op_list = op_list,
            op_l_vec = op_l_vec,
            αt_list = αt_list,
            δα_list = δα_list,
            dsf_cache = similar(ψ0.data, length(ψ0.data)^2),
            expv_cache = expv_cache,
            dsf_params = dsf_params,
            dsf_displace_cache_full = dsf_displace_cache_full,
        ),
    )

    cb_dsf = DiscreteCallback(_DSF_mesolve_Condition, _DSF_mesolve_Affect!, save_positions = (false, false))
    kwargs2 = (; kwargs...)
    kwargs2 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(cb_dsf, kwargs2.callback),)) :
        merge(kwargs2, (callback = cb_dsf,))

    return mesolveProblem(H₀, ψ0, tlist, c_ops₀; e_ops = e_ops₀, params = params2, kwargs2...)
end

@doc raw"""
    dsf_mesolve(H::Function,
        ψ0::QuantumObject,
        tlist::AbstractVector, c_ops::Function,
        op_list::Vector{TOl},
        α0_l::Vector{<:Number}=zeros(length(op_list)),
        dsf_params::NamedTuple=NamedTuple();
        alg::AbstractODEAlgorithm=DP5(),
        e_ops::Function=(op_list,p) -> Vector{TOl}([]),
        params::NamedTuple=NamedTuple(),
        δα_list::Vector{<:Number}=fill(0.2, length(op_list)),
        krylov_dim::Int=max(6, min(10, cld(length(ket2dm(ψ0).data), 4))),
        kwargs...)

Time evolution of an open quantum system using master equation and the Dynamical Shifted Fock algorithm.

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function dsf_mesolve(
        H::Function,
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        c_ops::Function,
        op_list::Union{AbstractVector, Tuple},
        α0_l::Vector{<:Number} = zeros(length(op_list)),
        dsf_params::NamedTuple = NamedTuple();
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Function = (op_list, p) -> (),
        params::NamedTuple = NamedTuple(),
        δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
        krylov_dim::Int = max(6, min(10, cld(length(ket2dm(ψ0).data), 4))),
        kwargs...,
    ) where {StateOpType <: Union{Ket, Operator}}
    dsf_prob = dsf_mesolveProblem(
        H,
        ψ0,
        tlist,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        e_ops = e_ops,
        params = params,
        δα_list = δα_list,
        krylov_dim = krylov_dim,
        kwargs...,
    )

    return mesolve(dsf_prob, alg)
end

function dsf_mesolve(
        H::Function,
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        op_list::Union{AbstractVector, Tuple},
        α0_l::Vector{<:Number} = zeros(length(op_list)),
        dsf_params::NamedTuple = NamedTuple();
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Function = (op_list, p) -> (),
        params::NamedTuple = NamedTuple(),
        δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
        krylov_dim::Int = max(6, min(10, cld(length(ket2dm(ψ0).data), 4))),
        kwargs...,
    ) where {StateOpType <: Union{Ket, Operator}}
    c_ops = op_list -> ()
    return dsf_mesolve(
        H,
        ψ0,
        tlist,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        alg = alg,
        e_ops = e_ops,
        params = params,
        δα_list = δα_list,
        krylov_dim = krylov_dim,
        kwargs...,
    )
end

# Dynamical Shifted Fock mcsolve

function _DSF_mcsolve_Condition(u, t, integrator)
    params = integrator.p
    op_list = params.op_list
    δα_list = params.δα_list

    ψt = u

    condition = false
    @inbounds for i in eachindex(op_list)
        op = op_list[i]
        δα = δα_list[i]
        Δα = dot(ψt, op.data, ψt) / dot(ψt, ψt)
        if δα < abs(Δα)
            condition = true
        end
    end
    return condition
end

function _DSF_mcsolve_Affect!(integrator)
    params = integrator.p
    op_list = params.op_list
    αt_list = params.αt_list
    δα_list = params.δα_list
    H = params.H_fun
    c_ops = params.c_ops_fun
    e_ops = params.e_ops_fun
    ψt = params.dsf_cache1
    dsf_cache = params.dsf_cache2
    expv_cache = params.expv_cache
    dsf_params = params.dsf_params
    dsf_displace_cache_full = params.dsf_displace_cache_full

    # e_ops0 = params.e_ops
    # c_ops0 = params.c_ops

    e_ops0 = _get_e_ops(integrator, SaveFuncMCSolve)
    c_ops0, c_ops0_herm = _mcsolve_get_c_ops(integrator)

    copyto!(ψt, integrator.u)
    normalize!(ψt)

    op_l_length = length(op_list)

    for i in eachindex(op_list)
        op = op_list[i]
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

            dsf_displace_cache_full.ops[i].λ.val = conj(Δα)
            dsf_displace_cache_full.ops[i + op_l_length].λ.val = -Δα

            αt_list[i] += Δα
        else
            dsf_displace_cache_full.ops[i].λ.val = 0
            dsf_displace_cache_full.ops[i + op_l_length].λ.val = 0
        end
    end

    copyto!(dsf_cache, integrator.u)
    arnoldi!(expv_cache, dsf_displace_cache_full, dsf_cache)
    expv!(integrator.u, expv_cache, 1, dsf_cache)

    op_l2 = op_list .+ αt_list
    e_ops2 = e_ops(op_l2, dsf_params)
    c_ops2 = c_ops(op_l2, dsf_params)

    ## By copying the data, we are assuming that the variables are Vectors and not Tuple
    @. e_ops0 = get_data(e_ops2)
    @. c_ops0 = get_data(c_ops2)
    c_ops0_herm .= map(op -> op' * op, c_ops0)

    H_nh = convert(eltype(ψt), 0.5im) * sum(c_ops0_herm)
    # By doing this, we are assuming that the system is time-independent and f is a ScaledOperator
    # of the form -1im * (H - H_nh)
    copyto!(integrator.f.f.L, H(op_l2, dsf_params).data - H_nh)
    return u_modified!(integrator, true)
end

function _dsf_mcsolve_prob_func(prob, i, repeat)
    params = prob.p

    prm = merge(
        params,
        (
            αt_list = copy(params.αt_list),
            dsf_cache1 = similar(params.dsf_cache1),
            dsf_cache2 = similar(params.dsf_cache2),
            expv_cache = copy(params.expv_cache),
            dsf_displace_cache_full = deepcopy(params.dsf_displace_cache_full), # This brutally copies also the MatrixOperators, and it is inefficient.
        ),
    )

    f = deepcopy(prob.f.f)

    # We need to deepcopy the callbacks because they contain the c_ops and e_ops, which are modified in the affect function. They also contain all the cache variables needed for mcsolve.
    cb = deepcopy(prob.kwargs[:callback])

    return remake(prob, f = f, p = prm, callback = cb)
end

function dsf_mcsolveEnsembleProblem(
        H::Function,
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Function,
        op_list::Union{AbstractVector, Tuple},
        α0_l::Vector{<:Number} = zeros(length(op_list)),
        dsf_params::NamedTuple = NamedTuple();
        e_ops::Function = (op_list, p) -> (),
        params::NamedTuple = NamedTuple(),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        krylov_dim::Int = min(5, cld(length(ψ0.data), 3)),
        progress_bar::Union{Bool, Val} = Val(true),
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    op_list = deepcopy(op_list)
    H₀ = H(op_list .+ α0_l, dsf_params)
    c_ops₀ = c_ops(op_list .+ α0_l, dsf_params)
    e_ops₀ = e_ops(op_list .+ α0_l, dsf_params)

    T = eltype(ψ0)

    αt_list = convert(Vector{T}, α0_l)
    expv_cache = arnoldi(H₀.data, ψ0.data, krylov_dim)

    dsf_displace_cache = sum(op -> ScalarOperator(one(T)) * MatrixOperator(op.data), op_list)
    dsf_displace_cache_dag = sum(op -> ScalarOperator(one(T)) * MatrixOperator(sparse(op.data')), op_list)
    dsf_displace_cache_full = dsf_displace_cache + dsf_displace_cache_dag

    params2 = merge(
        params,
        (
            H_fun = H,
            c_ops_fun = c_ops,
            e_ops_fun = e_ops,
            op_list = op_list,
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
        tlist,
        c_ops₀;
        e_ops = e_ops₀,
        params = params2,
        ntraj = ntraj,
        ensemblealg = ensemblealg,
        jump_callback = jump_callback,
        prob_func = _dsf_mcsolve_prob_func,
        progress_bar = progress_bar,
        kwargs2...,
    )
end

@doc raw"""
    dsf_mcsolve(H::Function,
        ψ0::QuantumObject,
        tlist::AbstractVector, c_ops::Function,
        op_list::Vector{TOl},
        α0_l::Vector{<:Number}=zeros(length(op_list)),
        dsf_params::NamedTuple=NamedTuple();
        alg::AbstractODEAlgorithm=DP5(),
        e_ops::Function=(op_list,p) -> Vector{TOl}([]),
        params::NamedTuple=NamedTuple(),
        δα_list::Vector{<:Real}=fill(0.2, length(op_list)),
        ntraj::Int=500,
        ensemblealg::EnsembleAlgorithm=EnsembleThreads(),
        jump_callback::LindbladJumpCallbackType=ContinuousLindbladJumpCallback(),
        krylov_dim::Int=max(6, min(10, cld(length(ket2dm(ψ0).data), 4))),
        progress_bar::Union{Bool,Val} = Val(true)
        kwargs...)

Time evolution of a quantum system using the Monte Carlo wave function method and the Dynamical Shifted Fock algorithm.

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function dsf_mcsolve(
        H::Function,
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector,
        c_ops::Function,
        op_list::Union{AbstractVector, Tuple},
        α0_l::Vector{<:Number} = zeros(length(op_list)),
        dsf_params::NamedTuple = NamedTuple();
        alg::AbstractODEAlgorithm = DP5(),
        e_ops::Function = (op_list, p) -> (),
        params::NamedTuple = NamedTuple(),
        δα_list::Vector{<:Real} = fill(0.2, length(op_list)),
        ntraj::Int = 500,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        jump_callback::TJC = ContinuousLindbladJumpCallback(),
        krylov_dim::Int = min(5, cld(length(ψ0.data), 3)),
        progress_bar::Union{Bool, Val} = Val(true),
        kwargs...,
    ) where {TJC <: LindbladJumpCallbackType}
    ens_prob_mc = dsf_mcsolveEnsembleProblem(
        H,
        ψ0,
        tlist,
        c_ops,
        op_list,
        α0_l,
        dsf_params;
        alg = alg,
        e_ops = e_ops,
        params = params,
        ntraj = ntraj,
        ensemblealg = ensemblealg,
        δα_list = δα_list,
        jump_callback = jump_callback,
        krylov_dim = krylov_dim,
        progress_bar = progress_bar,
        kwargs...,
    )

    return mcsolve(ens_prob_mc, alg, ntraj, ensemblealg)
end
