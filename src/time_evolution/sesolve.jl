export sesolveProblem, sesolve

function _save_func_sesolve(integrator)
    internal_params = integrator.p
    progr = internal_params.progr

    if !internal_params.is_empty_e_ops
        e_ops = internal_params.e_ops
        expvals = internal_params.expvals

        ψ = integrator.u
        _expect = op -> dot(ψ, op, ψ)
        @. expvals[:, progr.counter[]+1] = _expect(e_ops)
    end
    next!(progr)
    return u_modified!(integrator, false)
end

function _generate_sesolve_kwargs_with_callback(t_l, kwargs)
    cb1 = PresetTimeCallback(t_l, _save_func_sesolve, save_positions = (false, false))
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(kwargs.callback, cb1),)) :
        merge(kwargs, (callback = cb1,))

    return kwargs2
end

function _generate_sesolve_kwargs(e_ops, progress_bar::Val{true}, t_l, kwargs)
    return _generate_sesolve_kwargs_with_callback(t_l, kwargs)
end

function _generate_sesolve_kwargs(e_ops, progress_bar::Val{false}, t_l, kwargs)
    if e_ops isa Nothing
        return kwargs
    end
    return _generate_sesolve_kwargs_with_callback(t_l, kwargs)
end

@doc raw"""
    sesolveProblem(H,
        ψ0,
        tlist;
        alg=Tsit5()
        e_ops = nothing,
        params=NamedTuple(),
        progress_bar=Val(true),
        kwargs...)

Generates the ODEProblem for the Schrödinger time evolution of a quantum system:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle = -i \hat{H} |\psi(t)\rangle
```

# Arguments

- `H::Union{QuantumObject,Tuple}`: The Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: The initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: The time list of the evolution.
- `alg::OrdinaryDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::NamedTuple`: The parameters of the system.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: The keyword arguments passed to the `ODEProblem` constructor.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is empty, the default value of `saveat=tlist` (saving the states corresponding to `tlist`), otherwise, `saveat=[tlist[end]]` (only save the final state). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The `ODEProblem` for the Schrödinger time evolution of the system.
"""
function sesolveProblem(
    H::Union{QuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    ϕ0 = sparse_to_dense(_CType(ψ0), get_data(ψ0)) # Convert it to dense vector with complex element type

    t_l = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    H_evo = QobjEvo(H, -1im) # pre-multiply by -i
    isoper(H_evo) || throw(ArgumentError("The Hamiltonian must be an Operator."))
    check_dims(H_evo, ψ0)
    U = H_evo.data

    progr = ProgressBar(length(t_l), enable = getVal(progress_bar))

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(t_l))
        e_ops_data = ()
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        e_ops_data = get_data.(e_ops)
        is_empty_e_ops = isempty(e_ops)
    end

    p = (
        e_ops = e_ops_data,
        expvals = expvals,
        progr = progr,
        Hdims = H_evo.dims,
        times = t_l,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = is_empty_e_ops ? t_l : [t_l[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, makeVal(progress_bar), t_l, kwargs2)

    tspan = (t_l[1], t_l[end])
    return ODEProblem{true,FullSpecialize}(U, ϕ0, tspan, p; kwargs3...)
end

@doc raw"""
    sesolve(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Time evolution of a closed quantum system using the Schrödinger equation:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle = -i \hat{H} |\psi(t)\rangle
```

# Arguments

- `H::QuantumObject`: The Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: The initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: List of times at which to save the state of the system.
- `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::Union{Nothing,AbstractVector,Tuple}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
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
function sesolve(
    H::Union{QuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2}
    prob = sesolveProblem(H, ψ0, tlist; e_ops = e_ops, params = params, progress_bar = progress_bar, kwargs...)

    return sesolve(prob, alg)
end

function sesolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob, alg)

    ψt = map(ϕ -> QuantumObject(ϕ, type = Ket, dims = sol.prob.p.Hdims), sol.u)

    return TimeEvolutionSol(
        sol.prob.p.times,
        ψt,
        sol.prob.p.expvals,
        sol.retcode,
        sol.alg,
        sol.prob.kwargs[:abstol],
        sol.prob.kwargs[:reltol],
    )
end
