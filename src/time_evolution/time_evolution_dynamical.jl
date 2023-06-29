### DYNAMICAL FOCK DIMENSION ###

function _reduce_dims(QO::AbstractArray{T}, dims::Vector{<:Integer}, sel::AbstractVector, reduce::AbstractVector) where {T}
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

    ρmat
end

function _increase_dims(QO::AbstractArray{T}, dims::Vector{<:Integer}, sel::AbstractVector, increase::AbstractVector) where {T}
    rd = dims
    nd = length(rd)
    rd_new = zero(rd)
    rd_new[sel] .= increase
    @. rd_new = rd + rd_new

    if nd == 1
        ρmat = similar(QO, rd_new[1], rd_new[1])
        selectdim(ρmat, 1, rd[1]+1:rd_new[1]) .= 0
        selectdim(ρmat, 2, rd[1]+1:rd_new[1]) .= 0
        copyto!(view(ρmat, 1:rd[1], 1:rd[1]), QO)
    else
        ρmat2 = similar(QO, reverse!(repeat(rd_new, 2))...)
        ρmat = reshape(QO, reverse!(repeat(rd, 2))...)
        for i in eachindex(sel)
            selectdim(ρmat2, nd-sel[i]+1, rd[sel[i]]+1:rd_new[sel[i]]) .= 0
            selectdim(ρmat2, 2*nd-sel[i]+1, rd[sel[i]]+1:rd_new[sel[i]]) .= 0
        end
        copyto!(view(ρmat2, reverse!(repeat([1:n for n in rd], 2))...), ρmat)
        ρmat = reshape(ρmat2, prod(rd_new), prod(rd_new))
    end

    ρmat
end

_dfd_set_pillow = dim -> min(max(round(Int, 0.02 * dim), 1), 20)

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
    copyto!(dfd_ρt_cache, integrator.u)
    
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
    any(increase_list) || any(reduce_list)
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
    increase_list .= false
    reduce_list .= false
    push!(dim_list_evo_times, integrator.t)
    push!(dim_list_evo, dim_list)

    e_ops2 = map(op -> mat2vec(get_data(op)'), e_ops(dim_list))
    L = liouvillian(H(dim_list), c_ops(dim_list)).data

    resize!(integrator, size(L, 1))
    integrator.u .= mat2vec(ρt)
    integrator.p = merge(internal_params, (L = L, e_ops = e_ops2, 
                            dfd_ρt_cache = similar(integrator.u)))
end

function dfd_mesolveProblem(H::Function, ψ0::QuantumObject{<:AbstractArray{T1},StateOpType},
    t_l::AbstractVector, c_ops::Function, maxdims::Vector{T2};
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::Union{Nothing, Function}=nothing, 
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    tol_list::Vector{<:Number}=Float64[],
    kwargs...) where {T1,T2<:Integer,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    length(ψ0.dims) != length(maxdims) && throw(DimensionMismatch("'dim_list' and 'maxdims' do not have the same dimension."))

    e_ops === nothing && (e_ops = dim_list -> Vector{Vector{T1}}([]))

    dim_list = deepcopy(ψ0.dims)
    H₀ = H(dim_list)
    c_ops₀ = c_ops(dim_list)
    e_ops₀ = e_ops(dim_list)

    length(tol_list) != length(dim_list) && (tol_list = [1e-8 for d in dim_list])
    dim_list_evo_times = [0.0]
    dim_list_evo = [dim_list]
    reduce_list = zeros(Bool, length(dim_list))
    increase_list = zeros(Bool, length(dim_list))
    pillow_list = _dfd_set_pillow.(dim_list)

    params2 = merge(params, Dict(:H_fun => H, :c_ops_fun => c_ops, :e_ops_fun => e_ops,
                        :dim_list => dim_list, :maxdims => maxdims, :tol_list => tol_list,
                        :reduce_list => reduce_list, :increase_list => increase_list, :pillow_list => pillow_list,
                        :dim_list_evo_times => dim_list_evo_times, :dim_list_evo => dim_list_evo,
                        :dfd_ρt_cache => similar(ψ0.data, length(ψ0.data)^2)))

    cb_dfd = DiscreteCallback(_DFDIncreaseReduceCondition, _DFDIncreaseReduceAffect!, save_positions=(false, false))
    kwargs2 = kwargs
    kwargs2 = merge(kwargs2, haskey(kwargs2, :callback) ? 
                    Dict(:callback => CallbackSet(cb_dfd, kwargs2[:callback])) : Dict(:callback => cb_dfd))

    mesolveProblem(H₀, ψ0, t_l, c_ops₀; e_ops=e_ops₀, alg=alg, H_t=H_t, progress=progress,
                                params=params2, kwargs2...)
end

"""
    function dfd_mesolve(H::Function, ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function, maxdims::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
        e_ops::Union{Nothing, Function}=nothing, 
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        tol_list::Vector{<:Number}=Float64[],
        kwargs...)

Time evolution of an open quantum system using master equation, dynamically changing the dimension of the Hilbert subspaces.
"""
function dfd_mesolve(H::Function, ψ0::QuantumObject{<:AbstractArray{T1},StateOpType},
    t_l::AbstractVector, c_ops::Function, maxdims::Vector{T2};
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::Union{Nothing, Function}=nothing, 
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    tol_list::Vector{<:Number}=Float64[],
    kwargs...) where {T1,T2<:Integer,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    dfd_prob = dfd_mesolveProblem(H, ψ0, t_l, c_ops, maxdims; alg=alg, e_ops=e_ops, H_t=H_t, params=params,
                                    progress=progress, tol_list=tol_list, kwargs...)

    sol = solve(dfd_prob, alg)

    ρt = map(i -> QuantumObject(vec2mat(sol.u[i]), dims=sol.prob.p.dim_list_evo[searchsortedlast(sol.prob.p.dim_list_evo_times, sol.t[i])]), eachindex(sol.t))

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
        Δα = dot(op_vec, integrator.u)
        if δα < abs(Δα)
            condition = true
        end
    end
    condition
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
    L = internal_params.L
    dsf_cache = internal_params.dsf_cache

    for i in eachindex(op_l)
        op = op_l[i]
        op_vec = op_l_vec[i]
        αt = αt_list[i]
        δα = δα_list[i]
        Δα = dot(op_vec, integrator.u)

        if δα < abs(Δα)
            Dᵢ = exp(Δα*op' - conj(Δα)*op)
            copyto!(dsf_cache, integrator.u)
            mul!(integrator.u, sprepost(Dᵢ', Dᵢ).data, dsf_cache)

            αt_list[i] += Δα
        end
    end

    op_l2 = op_l .+ αt_list
    e_ops2 = e_ops(op_l2)
    _mat2vec_data = op -> mat2vec(get_data(op)')
    @. e_ops_vec = _mat2vec_data(e_ops2)
    copyto!(internal_params.L, liouvillian(H(op_l2), c_ops(op_l2)).data)
end

function dsf_mesolveProblem(H::Function,
    ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
    t_l::AbstractVector, c_ops::Function,
    op_list::AbstractVector,
    α0_l::Vector{<:Number}=zeros(length(op_list));
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::Union{Nothing, Function}=nothing,
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    δα_list::Vector{<:Real}=Float64[],
    kwargs...) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    e_ops === nothing && (e_ops = op_list -> Vector{QuantumObject{Matrix{ComplexF64}, OperatorQuantumObject}}([]))

    op_l = deepcopy(op_list)
    H₀ = H(op_l)
    c_ops₀ = c_ops(op_l)
    e_ops₀ = e_ops(op_l .+ α0_l)

    αt_list  = convert(Vector{T}, α0_l)
    length(δα_list) != length(op_l) ? δα_list = [0.3 for op in op_l] : nothing
    op_l_vec = map(op -> mat2vec(get_data(op)'), op_l)

    params2 = merge(params, Dict(:H_fun => H, :c_ops_fun => c_ops, :e_ops_fun => e_ops,
                    :op_l => op_l, :op_l_vec => op_l_vec, :αt_list => αt_list, :δα_list => δα_list,
                    :dsf_cache => similar(ψ0.data, length(ψ0.data)^2)))

    cb_dsf = DiscreteCallback(_DSF_mesolve_Condition, _DSF_mesolve_Affect!, save_positions=(false, false))
    kwargs2 = kwargs
    kwargs2 = merge(kwargs2, haskey(kwargs2, :callback) ? 
                    Dict(:callback => CallbackSet(cb_dsf, kwargs2[:callback])) : Dict(:callback => cb_dsf))

    mesolveProblem(H₀, ψ0, t_l, c_ops₀; e_ops=e_ops₀, alg=alg, H_t=H_t, progress=progress,
                                params=params2, kwargs2...)
end

"""
    function dsf_mesolve(H::Function,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function,
        op_list::AbstractVector,
        α0_l::Vector{<:Number}=zeros(length(op_list));
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
        e_ops::Union{Nothing, Function}=nothing,
        H_t::Union{Nothing,Function}=nothing,
        params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        progress::Bool=true,
        δα_list::Vector{<:Number}=Float64[],
        kwargs...)

Time evolution of an open quantum system using master equation and the Dynamical Shifted Fock algorithm.
"""
function dsf_mesolve(H::Function,
    ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
    t_l::AbstractVector, c_ops::Function,
    op_list::AbstractVector,
    α0_l::Vector{<:Number}=zeros(length(op_list));
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Vern7(),
    e_ops::Union{Nothing, Function}=nothing,
    H_t::Union{Nothing,Function}=nothing,
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    progress::Bool=true,
    δα_list::Vector{<:Real}=Float64[],
    kwargs...) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    dsf_prob = dsf_mesolveProblem(H, ψ0, t_l, c_ops, op_list, α0_l; alg=alg, e_ops=e_ops, H_t=H_t,
                                    params=params, progress=progress, δα_list=δα_list, kwargs...)
    
    return mesolve(dsf_prob; alg=alg, kwargs...)
end

# Dynamical Shifted Fock mcsolve

function _save_func_mcsolve_dsf(u, t, integrator)
    internal_params = integrator.p[1]
    save_it = internal_params["save_it"]
    op_l = internal_params["op_l"]
    αt_list = internal_params["αt_list"]
    e_ops = internal_params["e_ops"]
    expvals = internal_params["expvals"]
    ψ = normalize(u)
    expvals[:, save_it[]+1] .= map(op -> dot(ψ, op.data, ψ), e_ops(op_l .+ αt_list))
    save_it[]+=1
end

function LindbladJumpAffect_dsf!(integrator)
    ψ = integrator.u
    internal_params = integrator.p[1]
    op_l = internal_params["op_l"]
    αt_list = internal_params["αt_list"]
    c_ops = map(op -> op.data, internal_params["c_ops"](op_l .+ αt_list))

    if length(c_ops) == 1
        integrator.u = normalize!(c_ops[1] * ψ)
    else
        collaps_idx = 1
        r2 = rand()
        dp = 0
        c_op_ψ_l = Vector{Float64}(undef, length(c_ops))
        @inbounds for i in eachindex(c_ops)
            c_op_ψ = c_ops[i] * ψ
            res = real(dot(c_op_ψ, c_op_ψ))
            c_op_ψ_l[i] = res
            dp += res
        end
        prob = 0
        @inbounds for i in eachindex(c_ops)
            res = c_op_ψ_l[i]
            prob += res / dp
            if prob >= r2
                collaps_idx = i
                break
            end
        end
        integrator.u = normalize!(c_ops[collaps_idx] * ψ)
    end
    integrator.p[1]["random_n"] = rand()
end

function _DSF_mcsolve_Condition(u, t, integrator)
    internal_params = integrator.p[1]
    op_l = internal_params["op_l"]
    δα_list = internal_params["δα_list"]

    ψt = normalize(integrator.u)

    condition = false
    @inbounds for i in eachindex(op_l)
        op = op_l[i]
        δα = δα_list[i]
        Δα = dot(ψt, op.data, ψt)
        if δα < abs(Δα)
            condition = true
        end
    end
    condition
end

function _DSF_mcsolve_Affect!(integrator)
    internal_params = integrator.p[1]
    op_l = internal_params["op_l"]
    αt_list = internal_params["αt_list"]
    δα_list = internal_params["δα_list"]
    H = internal_params["H_fun"]
    H_eff = internal_params["H"]
    c_ops = internal_params["c_ops"]
    op1 = op_l[1]

    ψt = normalize(integrator.u)

    U = QuantumObject(spdiagm(ones(ComplexF64, size(op1, 1))), OperatorQuantumObject, op1.dims)
    @inbounds for i in eachindex(op_l)
        op = op_l[i]
        αt = αt_list[i]
        δα = δα_list[i]
        Δα = dot(ψt, op.data, ψt)
        
        if δα < abs(Δα)
            U *= exp(Δα*op' - conj(Δα)*op)
            αt_list[i] += Δα
        end
    end

    op_l2 = op_l .+ αt_list
    H_eff0 = H(op_l2).data
    for c_op in c_ops(op_l2)
        H_eff0 += -0.5im * c_op.data' * c_op.data
    end
    H_eff0 = -1im * H_eff0
    copyto!(H_eff, H_eff0)
    integrator.u = U.data' * integrator.u
end



function LindbladJumpAffect2!(integrator)
    ψ = integrator.u
    internal_params = integrator.p[1]
    c_ops = internal_params["c_ops"]

    if length(c_ops) == 1
        integrator.u = normalize!(c_ops[1] * ψ)
    else
        collaps_idx = 1
        r2 = rand()
        dp = 0
        c_op_ψ_l = Vector{Float64}(undef, length(c_ops))
        @inbounds for i in eachindex(c_ops)
            c_op_ψ = c_ops[i] * ψ
            res = real(dot(c_op_ψ, c_op_ψ))
            c_op_ψ_l[i] = res
            dp += res
        end
        prob = 0
        @inbounds for i in eachindex(c_ops)
            res = c_op_ψ_l[i]
            prob += res / dp
            if prob >= r2
                collaps_idx = i
                break
            end
        end
        integrator.u = normalize!(c_ops[collaps_idx] * ψ)
    end
    integrator.p[1]["random_n"] = rand()
end

function ContinuousLindbladJumpCallback2(interp_points::Int=0; affect_func::Function=LindbladJumpAffect2!)
    LindbladJumpCondition(u, t, integrator) = integrator.p[1]["random_n"] - real(dot(u, u))

    ContinuousCallback(LindbladJumpCondition, affect_func, nothing, interp_points=interp_points, save_positions=(false, false))
end

function DiscreteLindbladJumpCallback2(;affect_func::Function=LindbladJumpAffect2!)
    LindbladJumpCondition(u, t, integrator) = real(dot(u, u)) < integrator.p[1]["random_n"]

    DiscreteCallback(LindbladJumpCondition, affect_func, save_positions=(false, false))
end



"""
    dsf_mcsolve(H::Function, α0_l::Vector{<:Number},
        δ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
        t_l::AbstractVector, c_ops::Function, e_ops::Function, op_list::AbstractVector;
        δα_list::AbstractVector = [],
        n_traj::Int=1,
        batch_size::Int=min(Threads.nthreads(), n_traj),
        alg=AutoVern7(KenCarp4(autodiff=false)),
        ensemble_method=EnsembleThreads(),
        H_t=nothing,
        params::AbstractVector=[],
        progress::Bool=true,
        jump_interp_pts::Int=10,
        callbacks=[],
        kwargs...)

Time evolution of an open quantum system using quantum trajectories and the Dynamical Shifted Fock algorithm.
"""
function dsf_mcsolve(H::Function, α0_l::Vector{<:Number},
    δ0::QuantumObject{<:AbstractArray{T},KetQuantumObject},
    t_l::AbstractVector, c_ops::Function, e_ops::Function, op_list::AbstractVector;
    δα_list::AbstractVector = [],
    n_traj::Int=1,
    batch_size::Int=min(Threads.nthreads(), n_traj),
    alg=AutoVern7(KenCarp4(autodiff=false)),
    ensemble_method=EnsembleThreads(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    jump_interp_pts::Int=10,
    callbacks=[],
    kwargs...) where {T}

    op_l = map(i -> op_list[i] + α0_l[i], eachindex(α0_l))
    H(op_l).dims != δ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H(op_l).dims

    tspan = (t_l[1], t_l[end])
    e_ops_len = length(e_ops(op_l))

    δ0 = δ0.data

    progr = Progress(n_traj, showspeed=true, enabled=progress)
    channel = RemoteChannel(() -> Channel{Bool}(), 1)
    @async while take!(channel)
        next!(progr)
    end

    function prob_func(prob, i, repeat)
        op_l = deepcopy(op_list)
        αt_list  = convert(Vector{ComplexF64}, α0_l)
        length(δα_list) == 0 ? δα_list = [0.1 for op in op_l] : nothing
        H_eff = H(op_l .+ αt_list)
        for c_op in c_ops(op_l .+ αt_list)
            H_eff += -0.5im * c_op' * c_op
        end
        H_eff = -1im * H_eff.data
        expvals = Array{ComplexF64}(undef, e_ops_len, length(t_l))
        
        p = [Dict("H" => H_eff, "c_ops" => c_ops, "e_ops" => e_ops, "random_n" => rand(),
        "expvals" => expvals, "save_it" => Ref{Int32}(0),
        "op_l" => op_l, "αt_list" => αt_list, "δα_list" => δα_list,
        "H_fun" => H), params...]
        remake(prob, p=p)
    end
    function output_func(sol, i)
        put!(channel, true)
        if e_ops_len == 0
            res = map(ϕ -> QuantumObject(ϕ, dims=Hdims), sol.u)
        else
            res = sol.prob.p[1]["expvals"]
        end
        (res, false)
    end
    function reduction(u, batch, I)
        if e_ops_len == 0
            tmp = hcat(batch...)
            length(u) == 0 && return tmp, false
            res = hcat(u, tmp)
        else
            tmp = sum(cat(batch..., dims=3), dims=3)
            length(u) == 0 && return tmp, false
            res = sum(cat(u, tmp, dims=3), dims=3)
        end
        return res, false
    end

    is_time_dependent = !(H_t === nothing)
    if is_time_dependent
        dudt! = (du, u, p, t) -> mul!(du, p[1]["H"] - 1im * H_t(t).data, u)
    else
        dudt! = (du, u, p, t) -> mul!(du, p[1]["H"], u)
    end

    cb1 = FunctionCallingCallback(_save_func_mcsolve_dsf, funcat=t_l)
    cb2 = DiscreteCallback(_DSF_mcsolve_Condition, _DSF_mcsolve_Affect!, save_positions=(false,false))
    cb3 = jump_interp_pts == -1 ? DiscreteLindbladJumpCallback2(affect_func=LindbladJumpAffect_dsf!) : ContinuousLindbladJumpCallback2(jump_interp_pts, affect_func=LindbladJumpAffect_dsf!)
    cb4 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, cb3, cb4, callbacks...)

    prob = ODEProblem(dudt!, δ0, tspan, callback=cb; kwargs...)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, reduction=reduction)
    sol = solve(ensemble_prob, alg, ensemble_method, trajectories=n_traj, batch_size=batch_size)

    put!(channel, false)

    e_ops_len == 0 && return TimeEvolutionSol(Vector{Float64}(t_l), sol.u, [])

    e_ops_expect = dropdims(sum(sol.u, dims=3), dims=3) ./ n_traj

    return TimeEvolutionSol(Vector{Float64}(t_l), [], e_ops_expect)
end