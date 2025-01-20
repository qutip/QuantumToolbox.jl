export mesolveProblem, mesolve

_mesolve_make_L_QobjEvo(H::QuantumObject, c_ops) = QobjEvo(liouvillian(H, c_ops); type = SuperOperator)
_mesolve_make_L_QobjEvo(H::Union{QuantumObjectEvolution,Tuple}, c_ops) = liouvillian(QobjEvo(H), c_ops)

@doc raw"""
    mesolveProblem(
        H::Union{AbstractQuantumObject{HOpType},Tuple},
        ψ0::QuantumObject{StateOpType},
        tlist,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        inplace::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the ODEProblem for the master equation time evolution of an open quantum system:

```math
\frac{\partial \hat{\rho}(t)}{\partial t} = -i[\hat{H}, \hat{\rho}(t)] + \sum_n \mathcal{D}(\hat{C}_n) [\hat{\rho}(t)]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\hat{\rho}(t)] = \hat{C}_n \hat{\rho}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}(t) - \frac{1}{2} \hat{\rho}(t) \hat{C}_n^\dagger \hat{C}_n
```

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use the inplace version of the ODEProblem. The default is `Val(true)`. It is recommended to use `Val(true)` for better performance, but it is sometimes necessary to use `Val(false)`, for example when performing automatic differentiation using [Zygote.jl](https://github.com/FluxML/Zygote.jl).
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob::ODEProblem`: The ODEProblem for the master equation time evolution.
"""
function mesolveProblem(
    H::Union{AbstractQuantumObject{HOpType},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    tlist = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    L_evo = _mesolve_make_L_QobjEvo(H, c_ops)
    check_dimensions(L_evo, ψ0)

    T = Base.promote_eltype(L_evo, ψ0)

    ρ0 = sparse_to_dense(_CType(T), mat2vec(ket2dm(ψ0).data)) # Convert it to dense vector with complex element type
    L = L_evo.data

    is_empty_e_ops = (e_ops isa Nothing) ? true : isempty(e_ops)

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_se_me_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2, SaveFuncMESolve)

    tspan = (tlist[1], tlist[end])
    prob = ODEProblem{getVal(inplace),FullSpecialize}(L, ρ0, tspan, params; kwargs3...)

    return TimeEvolutionProblem(prob, tlist, L_evo.dimensions)
end

@doc raw"""
    mesolve(
        H::Union{AbstractQuantumObject{HOpType},Tuple},
        ψ0::QuantumObject{StateOpType},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        inplace::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Time evolution of an open quantum system using Lindblad master equation:

```math
\frac{\partial \hat{\rho}(t)}{\partial t} = -i[\hat{H}, \hat{\rho}(t)] + \sum_n \mathcal{D}(\hat{C}_n) [\hat{\rho}(t)]
```

where 

```math
\mathcal{D}(\hat{C}_n) [\hat{\rho}(t)] = \hat{C}_n \hat{\rho}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}(t) - \frac{1}{2} \hat{\rho}(t) \hat{C}_n^\dagger \hat{C}_n
```

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm for the ODE solver. The default value is `Tsit5()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use the inplace version of the ODEProblem. The default is `Val(true)`. It is recommended to use `Val(true)` for better performance, but it is sometimes necessary to use `Val(false)`, for example when performing automatic differentiation using [Zygote.jl](https://github.com/FluxML/Zygote.jl).
- `kwargs`: The keyword arguments for the ODEProblem.

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
    H::Union{AbstractQuantumObject{HOpType},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
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
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )

    return mesolve(prob, alg)
end

function mesolve(prob::TimeEvolutionProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob.prob, alg)

    ρt = map(ϕ -> QuantumObject(vec2mat(ϕ), type = Operator, dims = prob.dimensions), sol.u)

    return TimeEvolutionSol(
        prob.times,
        ρt,
        _se_me_sse_get_expvals(sol),
        sol.retcode,
        sol.alg,
        NamedTuple(sol.prob.kwargs).abstol,
        NamedTuple(sol.prob.kwargs).reltol,
    )
end
