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

sesolve_ti_dudt!(du, u, p, t) = mul!(du, p.U, u)
function sesolve_td_dudt!(du, u, p, t)
    mul!(du, p.U, u)
    H_t = p.H_t(t, p)
    return mul!(du, H_t, u, -1im, 1)
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
    sesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector;
        alg::OrdinaryDiffEqAlgorithm=Tsit5()
        e_ops::Union{Nothing,AbstractVector} = nothing,
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Union{Val,Bool}=Val(true),
        kwargs...)

Generates the ODEProblem for the Schrödinger time evolution of a quantum system:

```math
\frac{\partial}{\partial t} |\psi(t)\rangle = -i \hat{H} |\psi(t)\rangle
```

# Arguments

- `H::QuantumObject`: The Hamiltonian of the system ``\hat{H}``.
- `ψ0::QuantumObject`: The initial state of the system ``|\psi(0)\rangle``.
- `tlist::AbstractVector`: The time list of the evolution.
- `alg::OrdinaryDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::Union{Nothing,AbstractVector}`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::NamedTuple`: The parameters of the system.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: The keyword arguments passed to the `ODEProblem` constructor.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[tlist[end]]` (only save the final state), otherwise, `saveat=tlist` (saving the states corresponding to `tlist`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The `ODEProblem` for the Schrödinger time evolution of the system.
"""
function sesolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractVector{T2},KetQuantumObject},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2}
    H.dims != ψ0.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))

    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    is_time_dependent = !(H_t isa Nothing)
    progress_bar_val = makeVal(progress_bar)

    ϕ0 = _convert_u0(get_data(ψ0))

    t_l = convert(Vector{real(eltype(ϕ0))}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    U = -1im * get_data(H)
    progr = ProgressBar(length(t_l), enable = getVal(progress_bar_val))

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(t_l))
        e_ops2 = MT1[]
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
        e_ops2 = get_data.(e_ops)
        is_empty_e_ops = isempty(e_ops)
    end

    p = (
        U = U,
        e_ops = e_ops2,
        expvals = expvals,
        progr = progr,
        Hdims = H.dims,
        H_t = H_t,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = e_ops isa Nothing ? t_l : [t_l[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, progress_bar_val, t_l, kwargs2)

    dudt! = is_time_dependent ? sesolve_td_dudt! : sesolve_ti_dudt!

    tspan = (t_l[1], t_l[end])
    return ODEProblem{true,FullSpecialize}(dudt!, ϕ0, tspan, p; kwargs3...)
end

@doc raw"""
    sesolve(H::QuantumObject,
        ψ0::QuantumObject,
        tlist::AbstractVector;
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::Union{Nothing,AbstractVector} = nothing,
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
- `e_ops::Union{Nothing,AbstractVector}`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[tlist[end]]` (only save the final state), otherwise, `saveat=tlist` (saving the states corresponding to `tlist`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `sol::TimeEvolutionSol`: The solution of the time evolution. See also [`TimeEvolutionSol`](@ref)
"""
function sesolve(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractVector{T2},KetQuantumObject},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector} = nothing,
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {MT1<:AbstractMatrix,T2}
    prob = sesolveProblem(
        H,
        ψ0,
        tlist;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        progress_bar = makeVal(progress_bar),
        kwargs...,
    )

    return sesolve(prob, alg)
end

function sesolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob, alg)

    ψt = map(ϕ -> QuantumObject(ϕ, type = Ket, dims = sol.prob.p.Hdims), sol.u)

    return TimeEvolutionSol(
        sol.t,
        ψt,
        sol.prob.p.expvals,
        sol.retcode,
        sol.alg,
        sol.prob.kwargs[:abstol],
        sol.prob.kwargs[:reltol],
    )
end
