export mesolveProblem, mesolve

_mesolve_make_L_QobjEvo(H::Union{QuantumObject,Nothing}, c_ops) = QobjEvo(liouvillian(H, c_ops); type = SuperOperator())
_mesolve_make_L_QobjEvo(H::Union{QuantumObjectEvolution,Tuple}, c_ops) = liouvillian(QobjEvo(H), c_ops)
_mesolve_make_L_QobjEvo(H::Nothing, c_ops::Nothing) = throw(ArgumentError("Both H and
c_ops are Nothing. You are probably running the wrong function."))

@doc raw"""
    mesolveProblem(
        H::Union{AbstractQuantumObject,Tuple},
        ψ0::QuantumObject,
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
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use the inplace version of the ODEProblem. The default is `Val(true)`. It is recommended to use `Val(true)` for better performance, but it is sometimes necessary to use `Val(false)`, for example when performing automatic differentiation using [Zygote.jl](https://github.com/FluxML/Zygote.jl).
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- If `H` is an [`Operator`](@ref), `ψ0` is a [`Ket`](@ref) and `c_ops` is `Nothing`, the function will call [`sesolveProblem`](@ref) instead.
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
) where {HOpType<:Union{Operator,SuperOperator},StateOpType<:Union{Ket,Operator,OperatorKet}}
    (isoper(H) && isket(ψ0) && isnothing(c_ops)) && return sesolveProblem(
        H,
        ψ0,
        tlist;
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    tlist = _check_tlist(tlist, _float_type(ψ0))

    L_evo = _mesolve_make_L_QobjEvo(H, c_ops)
    check_dimensions(L_evo, ψ0)

    T = Base.promote_eltype(L_evo, ψ0)
    ρ0 = if isoperket(ψ0) # Convert it to dense vector with complex element type
        to_dense(_complex_float_type(T), copy(ψ0.data))
    else
        to_dense(_complex_float_type(T), mat2vec(ket2dm(ψ0).data))
    end
    L = L_evo.data

    kwargs2 = _merge_saveat(tlist, e_ops, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_se_me_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2, SaveFuncMESolve)

    tspan = (tlist[1], tlist[end])

    prob = ODEProblem{getVal(inplace),FullSpecialize}(L, ρ0, tspan, params; kwargs3...)

    return TimeEvolutionProblem(prob, tlist, L_evo.dimensions, (isoperket = Val(isoperket(ψ0)),))
end

@doc raw"""
    mesolve(
        H::Union{AbstractQuantumObject,Tuple},
        ψ0::QuantumObject,
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
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
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
- If `H` is an [`Operator`](@ref), `ψ0` is a [`Ket`](@ref) and `c_ops` is `Nothing`, the function will call [`sesolve`](@ref) instead.
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
) where {HOpType<:Union{Operator,SuperOperator},StateOpType<:Union{Ket,Operator,OperatorKet}}
    (isoper(H) && isket(ψ0) && isnothing(c_ops)) && return sesolve(
        H,
        ψ0,
        tlist;
        alg = alg,
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )

    # Move sensealg argument to solve for Enzyme.jl support.
    # TODO: Remove it when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed.
    sensealg = get(kwargs, :sensealg, nothing)
    kwargs_filtered = isnothing(sensealg) ? kwargs : Base.structdiff((; kwargs...), (sensealg = sensealg,))

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
        kwargs_filtered...,
    )

    # TODO: Remove sensealg when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed
    if isnothing(sensealg)
        return mesolve(prob, alg)
    else
        return mesolve(prob, alg; sensealg = sensealg)
    end
end

function mesolve(prob::TimeEvolutionProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)

    # No type instabilities since `isoperket` is a Val, and so it is known at compile time
    if getVal(prob.kwargs.isoperket)
        ρt = map(ϕ -> QuantumObject(ϕ, type = OperatorKet(), dims = prob.dimensions), sol.u)
    else
        ρt = map(ϕ -> QuantumObject(vec2mat(ϕ), type = Operator(), dims = prob.dimensions), sol.u)
    end

    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionSol(
        prob.times,
        sol.t,
        ρt,
        _get_expvals(sol, SaveFuncMESolve),
        sol.retcode,
        sol.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end
