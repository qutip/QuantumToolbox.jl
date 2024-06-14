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

@doc raw"""
    sesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5()
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Bool=true,
        kwargs...)

Generates the ODEProblem for the Schrödinger time evolution of a quantum system.

# Arguments

- `H::QuantumObject`: The Hamiltonian of the system.
- `ψ0::QuantumObject`: The initial state of the system.
- `t_l::AbstractVector`: The time list of the evolution.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: The algorithm used for the time evolution.
- `e_ops::AbstractVector`: The list of operators to be evaluated during the evolution.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: The time-dependent Hamiltonian of the system. If `nothing`, the Hamiltonian is time-independent.
- `params::NamedTuple`: The parameters of the system.
- `progress_bar::Bool`: Whether to show the progress bar.
- `kwargs...`: The keyword arguments passed to the `ODEProblem` constructor.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[t_l[end]]` (only save the final state), otherwise, `saveat=t_l` (saving the states corresponding to `t_l`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-5` and `abstol=1e-7`.
- For more details about `alg` and extra `kwargs`, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns

- `prob`: The `ODEProblem` for the Schrödinger time evolution of the system.
"""
function sesolveProblem(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{MT2,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Bool = true,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,MT2<:AbstractMatrix}
    H.dims != ψ0.dims && throw(DimensionMismatch("The two quantum objects are not of the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ϕ0 = get_data(ψ0)
    U = -1im * get_data(H)

    progr = ProgressBar(length(t_l), enable = progress_bar)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = get_data.(e_ops)
    is_empty_e_ops = isempty(e_ops)

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

    saveat = is_empty_e_ops ? t_l : [t_l[end]]
    default_values = (abstol = 1e-7, reltol = 1e-5, saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    if !isempty(e_ops) || progress_bar
        cb1 = PresetTimeCallback(t_l, _save_func_sesolve, save_positions = (false, false))
        kwargs2 =
            haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(kwargs2.callback, cb1),)) :
            merge(kwargs2, (callback = cb1,))
    end

    tspan = (t_l[1], t_l[end])
    return _sesolveProblem(U, ϕ0, tspan, alg, Val(is_time_dependent), p; kwargs2...)
end

function _sesolveProblem(
    U::AbstractMatrix{<:T1},
    ϕ0::AbstractVector{<:T2},
    tspan::Tuple,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    is_time_dependent::Val{false},
    p;
    kwargs...,
) where {T1,T2}
    return ODEProblem{true,SciMLBase.FullSpecialize}(sesolve_ti_dudt!, ϕ0, tspan, p; kwargs...)
end

function _sesolveProblem(
    U::AbstractMatrix{<:T1},
    ϕ0::AbstractVector{<:T2},
    tspan::Tuple,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    is_time_dependent::Val{true},
    p;
    kwargs...,
) where {T1,T2}
    return ODEProblem{true,SciMLBase.FullSpecialize}(sesolve_td_dudt!, ϕ0, tspan, p; kwargs...)
end

@doc raw"""
    sesolve(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector;
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Bool=true,
        kwargs...)

Time evolution of a closed quantum system using the Schrödinger equation:

```math
\frac{\partial \psi(t)}{\partial t}=-i H \psi(t)
```

# Arguments

- `H::QuantumObject`: Hamiltonian of the system.
- `ψ0::QuantumObject`: Initial state of the system.
- `t_l::AbstractVector`: List of times at which to save the state of the system.
- `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
- `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
- `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
- `params::NamedTuple`: Dictionary of parameters to pass to the solver.
- `progress_bar::Bool`: Whether to show the progress bar.
- `kwargs...`: Additional keyword arguments to pass to the solver.

# Notes

- The states will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[t_l[end]]` (only save the final state), otherwise, `saveat=t_l` (saving the states corresponding to `t_l`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-5` and `abstol=1e-7`.
- For more details about `alg` and extra `kwargs`, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns

- `sol::TimeEvolutionSol`: The solution of the time evolution. See also [`TimeEvolutionSol`](@ref)
"""
function sesolve(
    H::QuantumObject{MT1,OperatorQuantumObject},
    ψ0::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
    t_l::AbstractVector;
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{MT2,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Bool = true,
    kwargs...,
) where {MT1<:AbstractMatrix,T2,MT2<:AbstractMatrix}
    prob = sesolveProblem(
        H,
        ψ0,
        t_l;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        progress_bar = progress_bar,
        kwargs...,
    )

    return sesolve(prob, alg)
end

function sesolve(prob::ODEProblem, alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob, alg)
    ψt =
        isempty(sol.prob.kwargs[:saveat]) ? QuantumObject[] :
        map(ϕ -> QuantumObject(ϕ, type = Ket, dims = sol.prob.p.Hdims), sol.u)

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
