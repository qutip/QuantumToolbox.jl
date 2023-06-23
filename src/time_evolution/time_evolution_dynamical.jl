### DYNAMICAL FOCK DIMENSION ###

function _save_func_dfd_mesolve(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params["progr"]
    e_ops = internal_params["e_ops"]
    dim_list = internal_params["dim_list"]
    expvals = internal_params["expvals"]
    
    ρt = reshape(copy(u), prod(dim_list), :)
    expvals[:, progr.counter+1] .= map(op -> tr(op.data * ρt), e_ops(dim_list))
    next!(progr)
end

function _reduce_dims(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::AbstractVector, reduce::AbstractVector) where {T,OpType<:OperatorQuantumObject}
    rd = QO.dims
    nd = length(rd)
    reduce_l = zero(rd)
    reduce_l[sel] .= reduce
    rd_new = rd .- reduce_l

    if nd == 1
        ρmat = similar(QO.data, repeat(rd_new, 2)...)
        copyto!(ρmat, view(QO.data, repeat([1:n for n in rd_new], 2)...))
    else
        ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
        ρmat2 = similar(QO.data, repeat(rd_new, 2)...)
        copyto!(ρmat2, view(ρmat, repeat([1:n for n in rd_new], 2)...))
        ρmat = reshape(PermutedDimsArray(ρmat2, length(size(ρmat2)):-1:1), prod(rd_new), prod(rd_new))
    end

    QuantumObject(ρmat, OperatorQuantumObject, rd_new)
end

function _increase_dims(QO::QuantumObject{<:AbstractArray{T},OpType}, sel::AbstractVector, increase::AbstractVector) where {T,OpType<:OperatorQuantumObject}
    rd = QO.dims
    nd = length(rd)
    incr_l = zero(rd)
    incr_l[sel] .= increase
    rd_new = rd .+ incr_l

    if nd == 1
        ρmat = similar(QO.data, repeat(rd_new, 2)...)
        selectdim(ρmat, sel[1], rd[1]+1:rd_new[1]) .= 0
        selectdim(ρmat, sel[1] + nd, rd[1]+1:rd_new[1]) .= 0
        copyto!(view(ρmat, repeat([1:n for n in rd], 2)...), QO.data)
    else
        ρmat = row_major_reshape(QO.data, repeat(rd, 2)...)
        ρmat2 = similar(QO.data, repeat(rd_new, 2)...)
        for i in eachindex(sel)
            selectdim(ρmat2, sel[i], rd[i]+1:rd_new[i]) .= 0
            selectdim(ρmat2, sel[i] + nd, rd[i]+1:rd_new[i]) .= 0
        end
        copyto!(view(ρmat2, repeat([1:n for n in rd], 2)...), ρmat)
        ρmat = reshape(PermutedDimsArray(ρmat2, length(size(ρmat2)):-1:1), prod(rd_new), prod(rd_new))
    end

    QuantumObject(ρmat, OperatorQuantumObject, rd_new)
end

function _DFDIncreaseReduceCondition(u, t, integrator)
    internal_params = integrator.p[1]
    dim_list = internal_params["dim_list"]
    maxdims = internal_params["maxdims"]
    tol_list = internal_params["tol_list"]
    increase_list = internal_params["increase_list"]
    reduce_list = internal_params["reduce_list"]
    pillow_list = internal_params["pillow_list"]
    # Here I use copy(u) otherwise I have an error.
    ρt = QuantumObject(reshape(copy(u), prod(dim_list), prod(dim_list)), OperatorQuantumObject, dim_list)

    condition = false
    @inbounds for i in 1:length(dim_list)
        maxdim_i = maxdims[i]
        dim_i = dim_list[i]
        if dim_i < maxdim_i && dim_i > 2 && maxdim_i != 0
            ρi = ptrace(ρt, [i]).data
            pillow_i = min(max(round(Int, 0.02 * dim_i), 1), 20)
            pillow_list[i] = pillow_i
            @views res = sum(abs.(ρi[diagind(ρi)[end-pillow_i:end]])) * sqrt(dim_i) / pillow_i
            if res > tol_list[i]
                push!(increase_list, i)
                condition = true
            elseif res < tol_list[i] * 1e-2 && dim_i > 3
                push!(reduce_list, i)
                condition = true
            end
        end
    end
    condition
end

function _DFDIncreaseReduceAffect!(integrator)
    internal_params = integrator.p[1]
    H = internal_params["H"]
    c_ops = internal_params["c_ops"]
    dim_list = internal_params["dim_list"]
    increase_list = internal_params["increase_list"]
    reduce_list = internal_params["reduce_list"]
    pillow_list = internal_params["pillow_list"]
    # Here I use copy(integrator.u) otherwise I have an error.
    ρt = QuantumObject(reshape(copy(integrator.u), prod(dim_list), prod(dim_list)), OperatorQuantumObject, dim_list)

    @views pillow_increase = pillow_list[increase_list]
    @views pillow_reduce = pillow_list[reduce_list]
    if length(increase_list) > 0
        ρt = _increase_dims(ρt, increase_list, pillow_increase)
        dim_list[increase_list] .+= pillow_increase
    end
    if length(reduce_list) > 0
        ρt = _reduce_dims(ρt, reduce_list, pillow_reduce)
        dim_list[reduce_list] .-= pillow_reduce
    end

    deleteat!(increase_list, 1:length(increase_list))
    deleteat!(reduce_list, 1:length(reduce_list))

    L = liouvillian(H(dim_list), c_ops(dim_list)).data
    newsize = size(L, 1)
    resize!(integrator, newsize)
    internal_params["L"] = L
    integrator.u = reshape(ρt.data, newsize)
end

"""
    dfd_mesolve(H::Function, ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::Function, e_ops::Function, maxdims::AbstractVector;
        tol_list::AbstractVector=[],
        alg = Vern7(),
        params::AbstractVector=[],
        progress::Bool=true,
        callbacks::AbstractVector=[],
        kwargs...)

Time evolution of an open quantum system using master equation, dynamically changing the dimension of the Hilbert subspaces.
"""
function dfd_mesolve(H::Function, ψ0::QuantumObject{<:AbstractArray{T},StateOpType},
    t_l::AbstractVector, c_ops::Function, e_ops::Function, maxdims::AbstractVector;
    tol_list::AbstractVector=[],
    alg=Vern7(),
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks::AbstractVector=[],
    kwargs...) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    dim_list = copy(ψ0.dims)
    length(dim_list) != length(maxdims) && throw(DimensionMismatch("'dim_list' and 'maxdims' do not have the same dimension."))
    H(dim_list).dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ρ0 = isket(ψ0) ? reshape(ψ0.data * ψ0.data', length(ψ0)^2) : reshape(ψ0.data, length(ψ0))

    L = liouvillian(H(dim_list), c_ops(dim_list)).data

    length(tol_list) == 0 && (tol_list = [1e-8 for d in dim_list])
    reduce_list = Int16[]
    increase_list = Int16[]
    pillow_list = ones(Int16, length(dim_list))

    cb1 = FunctionCallingCallback(_save_func_dfd_mesolve, funcat=t_l)
    cb2 = AutoAbstol(false; init_curmax=0.0)
    cb_increasereduce = DiscreteCallback(_DFDIncreaseReduceCondition, _DFDIncreaseReduceAffect!, save_positions=(false, false))
    cb = CallbackSet(cb1, cb2, cb_increasereduce, callbacks...)

    expvals = Array{ComplexF64}(undef, length(e_ops(dim_list)), length(t_l))
    p = [Dict("L" => L, "H" => H, "c_ops" => c_ops, "e_ops" => e_ops,
        "dim_list" => dim_list, "maxdims" => maxdims, "tol_list" => tol_list,
        "reduce_list" => reduce_list, "increase_list" => increase_list, "pillow_list" => pillow_list,
        "expvals" => expvals, "progr" => progr), params...]

    dudt! = (du, u, p, t) -> mul!(du, p[1]["L"], u)
    prob = ODEProblem(dudt!, ρ0, tspan, p; kwargs...)
    sol = solve(prob, alg, callback=cb)

    ρt = [QuantumObject(reshape(ϕ, isqrt(length(ϕ)), :)) for ϕ in sol.u]

    TimeEvolutionSol(sol.t, ρt, sol.prob.p[1]["expvals"])
end

# Dynamical Shifted Fock mesolve

function _save_func_mesolve_dsf(u, t, integrator)
    internal_params = integrator.p[1]
    progr = internal_params["progr"]
    op_l = internal_params["op_l"]
    αt_list = internal_params["αt_list"]
    e_ops = internal_params["e_ops"]
    expvals = internal_params["expvals"]
    op1 = op_l[1]
    
    expvals[:, progr.counter+1] .= map(op -> tr(op.data * reshape(u, isqrt(length(u)), :)), e_ops(op_l .+ αt_list))
    next!(progr)
end

function _DSF_mesolve_Condition(u, t, integrator)
    internal_params = integrator.p[1]
    op_l = internal_params["op_l"]
    δα_list = internal_params["δα_list"]
    op1 = op_l[1]
    # integrator.u or just u?
    ρt = reshape(integrator.u, size(op1)...)

    condition = false
    @inbounds for i in eachindex(op_l)
        op = op_l[i]
        δα = δα_list[i]
        Δα = tr(op.data * ρt)
        if δα < abs(Δα)
            condition = true
        end
    end
    condition
end

function _DSF_mesolve_Affect!(integrator)
    internal_params = integrator.p[1]
    op_l = internal_params["op_l"]
    αt_list = internal_params["αt_list"]
    δα_list = internal_params["δα_list"]
    H = internal_params["H"]
    c_ops = internal_params["c_ops"]
    L = internal_params["L"]
    op1 = op_l[1]

    ρt  = reshape(integrator.u, size(op1)...)

    U = QuantumObject(spdiagm(ones(ComplexF64, size(op1, 1))), OperatorQuantumObject, op1.dims)
    @inbounds for i in eachindex(op_l)
        op = op_l[i]
        αt = αt_list[i]
        δα = δα_list[i]
        Δα = tr(op.data * ρt)

        if δα < abs(Δα)
            U *= exp(Δα*op' - conj(Δα)*op)
            αt_list[i] += Δα
        end
    end

    op_l2 = op_l .+ αt_list
    copyto!(L, liouvillian(H(op_l2), c_ops(op_l2)).data)
    integrator.u = sprepost(U', U).data * integrator.u
end

"""
    dsf_mesolve(H::Function, α0_l::Vector{<:Number},
        ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
        t_l::AbstractVector, c_ops::Function, e_ops::Function, op_list::AbstractVector;
        δα_list::AbstractVector = [],
        alg = Vern7(),
        params::AbstractVector=[],
        progress = true,
        callbacks = [],
        kwargs...)

Time evolution of an open quantum system using master equation and the Dynamical Shifted Fock algorithm.
"""
function dsf_mesolve(H::Function, α0_l::Vector{<:Number},
    ψ0::QuantumObject{<:AbstractArray{T}, StateOpType},
    t_l::AbstractVector, c_ops::Function, e_ops::Function, op_list::AbstractVector;
    δα_list::AbstractVector = [],
    alg = Vern7(),
    params::AbstractVector=[],
    progress = true,
    callbacks = [],
    kwargs...) where {T,StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}

    op_l = deepcopy(op_list)
    H(op_l).dims != ψ0.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    Hdims = H(op_l).dims

    tspan = (t_l[1], t_l[end])

    progr = Progress(length(t_l), showspeed=true, enabled=progress)

    ρ0 = isket(ψ0) ? reshape(ψ0.data * ψ0.data', length(ψ0)^2) : reshape(ψ0.data, length(ψ0))

    L = liouvillian(H(op_l), c_ops(op_l)).data

    αt_list  = convert(Vector{ComplexF64}, α0_l)
    length(δα_list) == 0 ? δα_list = [0.1 for op in op_l] : nothing
    expvals = Array{ComplexF64}(undef, length(e_ops(op_l)), length(t_l))
    p = [Dict("L" => L, "H" => H, "c_ops" => c_ops, "e_ops" => e_ops,
    "op_l" => op_l, "αt_list" => αt_list, "δα_list" => δα_list,
    "expvals" => expvals, "progr" => progr), params...]

    cb1 = FunctionCallingCallback(_save_func_mesolve_dsf, funcat=t_l)
    cb2 = DiscreteCallback(_DSF_mesolve_Condition, _DSF_mesolve_Affect!, save_positions=(false,false))
    cb3 = AutoAbstol(false; init_curmax=0.0)
    cb = CallbackSet(cb1, cb2, cb3, callbacks...)

    dudt! = (du, u, p, t) -> mul!(du, p[1]["L"], u)
    prob = ODEProblem(dudt!, ρ0, tspan, p; kwargs...)
    sol = solve(prob, alg, callback=cb)

    ρt_len = isqrt(length(sol.u[1]))
    ρt_len == prod(Hdims) ? ρt = [QuantumObject(reshape(ϕ, ρt_len, ρt_len), dims=Hdims) for ϕ in sol.u] : ρt = sol.u

    TimeEvolutionSol(sol.t, ρt, sol.prob.p[1]["expvals"])
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