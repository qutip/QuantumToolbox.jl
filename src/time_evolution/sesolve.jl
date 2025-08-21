export sesolveProblem, sesolve

_sesolve_make_U_QobjEvo(
    H::QuantumObjectEvolution{Operator,DimsType,<:MatrixOperator},
) where {DimsType<:AbstractDimensions} =
    QobjEvo(MatrixOperator(-1im * H.data.A), dims = H.dimensions, type = Operator())
_sesolve_make_U_QobjEvo(H::QuantumObject) =
    QobjEvo(MatrixOperator(-1im * H.data), dims = H.dimensions, type = Operator())
_sesolve_make_U_QobjEvo(H::Union{QuantumObjectEvolution,Tuple}) = QobjEvo(H, -1im)

@doc raw"""
    sesolveProblem(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        inplace::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the ODEProblem for the Schrödinger time evolution of a quantum system:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle = -i \hat{H} |\psi(t)\rangle
```

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
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

- `prob`: The [`TimeEvolutionProblem`](@ref) containing the `ODEProblem` for the Schrödinger time evolution of the system.
"""
function sesolveProblem(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
)
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    tlist = _check_tlist(tlist, _float_type(ψ0))

    H_evo = _sesolve_make_U_QobjEvo(H) # Multiply by -i
    isoper(H_evo) || throw(ArgumentError("The Hamiltonian must be an Operator."))
    check_dimensions(H_evo, ψ0)

    T = Base.promote_eltype(H_evo, ψ0)
    ψ0 = to_dense(_complex_float_type(T), get_data(ψ0)) # Convert it to dense vector with complex element type
    U = H_evo.data

    kwargs2 = _merge_saveat(tlist, e_ops, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_se_me_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2, SaveFuncSESolve)

    tspan = (tlist[1], tlist[end])

    prob = ODEProblem{getVal(inplace),FullSpecialize}(U, ψ0, tspan, params; kwargs3...)

    return TimeEvolutionProblem(prob, tlist, H_evo.dimensions)
end

@doc raw"""
    sesolve(
        H::Union{AbstractQuantumObject{Operator},Tuple},
        ψ0::QuantumObject{Ket},
        tlist::AbstractVector;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        inplace::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Time evolution of a closed quantum system using the Schrödinger equation:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle = -i \hat{H} |\psi(t)\rangle
```

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state of the system ``|\psi(0)\rangle``.
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `alg`: The algorithm for the ODE solver. The default is `Tsit5()`.
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
function sesolve(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
)

    # Move sensealg argument to solve for Enzyme.jl support.
    # TODO: Remove it when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed.
    sensealg = get(kwargs, :sensealg, nothing)
    kwargs_filtered = isnothing(sensealg) ? kwargs : Base.structdiff((; kwargs...), (sensealg = sensealg,))

    prob = sesolveProblem(
        H,
        ψ0,
        tlist;
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs_filtered...,
    )

    # TODO: Remove it when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed.
    if isnothing(sensealg)
        return sesolve(prob, alg)
    else
        return sesolve(prob, alg; sensealg = sensealg)
    end
end

function sesolve(prob::TimeEvolutionProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)

    ψt = map(ϕ -> QuantumObject(ϕ, type = Ket(), dims = prob.dimensions), sol.u)

    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionSol(
        prob.times,
        sol.t,
        ψt,
        _get_expvals(sol, SaveFuncSESolve),
        sol.retcode,
        sol.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end
