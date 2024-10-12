export mesolveProblem, mesolve

function _save_func_mesolve(integrator)
    internal_params = integrator.p
    progr = internal_params.progr

    if !internal_params.is_empty_e_ops
        expvals = internal_params.expvals
        e_ops = internal_params.e_ops
        # This is equivalent to tr(op * ρ), when both are matrices.
        # The advantage of using this convention is that I don't need
        # to reshape u to make it a matrix, but I reshape the e_ops once.

        ρ = integrator.u
        _expect = op -> dot(op, ρ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    next!(progr)
    return u_modified!(integrator, false)
end

mesolve_ti_dudt!(du, u, p, t) = mul!(du, p.L, u)
function mesolve_td_dudt!(du, u, p, t)
    mul!(du, p.L, u)
    L_t = p.H_t(t, p)
    return mul!(du, L_t, u, 1, 1)
end

_generate_mesolve_e_op(op) = mat2vec(adjoint(get_data(op)))

function _generate_mesolve_kwargs_with_callback(t_l, kwargs)
    cb1 = PresetTimeCallback(t_l, _save_func_mesolve, save_positions = (false, false))
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(kwargs.callback, cb1),)) :
        merge(kwargs, (callback = cb1,))

    return kwargs2
end

function _generate_mesolve_kwargs(e_ops, progress_bar::Val{true}, t_l, kwargs)
    return _generate_mesolve_kwargs_with_callback(t_l, kwargs)
end

function _generate_mesolve_kwargs(e_ops, progress_bar::Val{false}, t_l, kwargs)
    if e_ops isa Nothing
        return kwargs
    end
    return _generate_mesolve_kwargs_with_callback(t_l, kwargs)
end

@doc raw"""
    mesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector, 
        c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple}=nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Generates the ODEProblem for the master equation time evolution of an open quantum system:

```math
\frac{\partial \hat{\rho}(t)}{\partial t} = -i[\hat{H}, \hat{\rho}(t)] + \sum_n \mathcal{D}(\hat{C}_n) [\hat{\rho}(t)]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\hat{\rho}(t)] = \hat{C}_n \hat{\rho}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}(t) - \frac{1}{2} \hat{\rho}(t) \hat{C}_n^\dagger \hat{C}_n
```

# Arguments

- `H::QuantumObject`: The Hamiltonian ``\hat{H}`` or the Liouvillian of the system.
- `ψ0::QuantumObject`: The initial state of the system.
- `tlist::AbstractVector`: The time list of the evolution.
- `c_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: The list of the collapse operators ``\{\hat{C}_n\}_n``.
- `alg::OrdinaryDiffEqAlgorithm=Tsit5()`: The algorithm used for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: The list of the operators for which the expectation values are calculated.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing`: The time-dependent Hamiltonian or Liouvillian.
- `params::NamedTuple=NamedTuple()`: The parameters of the time evolution.
- `progress_bar::Union{Val,Bool}=Val(true)`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::ODEProblem`: The ODEProblem for the master equation time evolution.
"""
function mesolveProblem(
    H::QuantumObject{MT1,HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    tlist,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
    MT1<:AbstractMatrix,
    T2,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    check_dims(H, ψ0)

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    is_time_dependent = !(H_t isa Nothing)

    ρ0 = sparse_to_dense(_CType(ψ0), mat2vec(ket2dm(ψ0).data)) # Convert it to dense vector with complex element type

    t_l = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    L = liouvillian(H, c_ops).data
    progr = ProgressBar(length(t_l), enable = getVal(progress_bar))

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(t_l))
        e_ops2 = mat2vec(MT1)[]
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        e_ops2 = [_generate_mesolve_e_op(op) for op in e_ops]
        is_empty_e_ops = isempty(e_ops)
    end

    if H_t isa TimeDependentOperatorSum
        H_t = liouvillian(H_t)
    end

    p = (
        L = L,
        progr = progr,
        Hdims = H.dims,
        e_ops = e_ops2,
        expvals = expvals,
        H_t = H_t,
        times = t_l,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = is_empty_e_ops ? t_l : [t_l[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_mesolve_kwargs(e_ops, makeVal(progress_bar), t_l, kwargs2)

    dudt! = is_time_dependent ? mesolve_td_dudt! : mesolve_ti_dudt!

    tspan = (t_l[1], t_l[end])
    return ODEProblem{true,FullSpecialize}(dudt!, ρ0, tspan, p; kwargs3...)
end

@doc raw"""
    mesolve(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector, 
        c_ops::Union{Nothing,AbstractVector,Tuple}=nothing;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple}=nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Time evolution of an open quantum system using Lindblad master equation:

```math
\frac{\partial \hat{\rho}(t)}{\partial t} = -i[\hat{H}, \hat{\rho}(t)] + \sum_n \mathcal{D}(\hat{C}_n) [\hat{\rho}(t)]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\hat{\rho}(t)] = \hat{C}_n \hat{\rho}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}(t) - \frac{1}{2} \hat{\rho}(t) \hat{C}_n^\dagger \hat{C}_n
```

# Arguments

- `H::QuantumObject`: The Hamiltonian ``\hat{H}`` or the Liouvillian of the system.
- `ψ0::QuantumObject`: The initial state of the system.
- `tlist::AbstractVector`: The time list of the evolution.
- `c_ops::Union{Nothing,AbstractVector,Tuple}=nothing`: The list of the collapse operators ``\{\hat{C}_n\}_n``.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Named Tuple of parameters to pass to the solver.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionSol`: The solution of the time evolution. See also [`TimeEvolutionSol`](@ref)
"""
function mesolve(
    H::QuantumObject{MT1,HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
    MT1<:AbstractMatrix,
    T2,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    prob = mesolveProblem(
        H,
        ψ0,
        tlist,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        progress_bar = progress_bar,
        kwargs...,
    )

    return mesolve(prob, alg)
end

function mesolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob, alg)

    ρt = map(ϕ -> QuantumObject(vec2mat(ϕ), type = Operator, dims = sol.prob.p.Hdims), sol.u)

    return TimeEvolutionSol(
        sol.prob.p.times,
        ρt,
        sol.prob.p.expvals,
        sol.retcode,
        sol.alg,
        sol.prob.kwargs[:abstol],
        sol.prob.kwargs[:reltol],
    )
end
