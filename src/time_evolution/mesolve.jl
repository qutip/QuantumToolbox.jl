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

_generate_mesolve_e_op(op) = mat2vec(adjoint(get_data(op)))

function _generate_mesolve_kwargs_with_callback(tlist, kwargs)
    cb1 = PresetTimeCallback(tlist, _save_func_mesolve, save_positions = (false, false))
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(kwargs.callback, cb1),)) :
        merge(kwargs, (callback = cb1,))

    return kwargs2
end

function _generate_mesolve_kwargs(e_ops, progress_bar::Val{true}, tlist, kwargs)
    return _generate_mesolve_kwargs_with_callback(tlist, kwargs)
end

function _generate_mesolve_kwargs(e_ops, progress_bar::Val{false}, tlist, kwargs)
    if e_ops isa Nothing
        return kwargs
    end
    return _generate_mesolve_kwargs_with_callback(tlist, kwargs)
end

_mesolve_make_L_QobjEvo(H::QuantumObject, c_ops) = QobjEvo(liouvillian(H, c_ops); type = SuperOperator)
function _mesolve_make_L_QobjEvo(H::Tuple, c_ops)
    c_ops isa Nothing && return QobjEvo(H)
    return QobjEvo((H..., mapreduce(op -> lindblad_dissipator(op), +, c_ops)); type = SuperOperator, f = liouvillian)
end
_mesolve_make_L_QobjEvo(H::QuantumObjectEvolution{DT,OperatorQuantumObject}, c_ops) where {DT<:AbstractSciMLOperator} =
    throw(
        ArgumentError(
            "This function does not support the data type of time-dependent Operator `H` currently. Try to provide `H` as a time-dependent SuperOperator or Tuple instead.",
        ),
    )
function _mesolve_make_L_QobjEvo(
    H::QuantumObjectEvolution{DT,SuperOperatorQuantumObject},
    c_ops,
) where {DT<:AbstractSciMLOperator}
    c_ops isa Nothing && return H
    return H + QobjEvo((mapreduce(op -> lindblad_dissipator(op), +, c_ops)))
end

@doc raw"""
    mesolveProblem(
        H::Union{AbstractQuantumObject{DT1,HOpType},Tuple},
        ψ0::QuantumObject{DT2,StateOpType},
        tlist,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        progress_bar::Union{Val,Bool} = Val(true),
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
- `params`: `NamedTuple` of parameters to pass to the solver.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
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
    H::Union{AbstractQuantumObject{DT1,HOpType},Tuple},
    ψ0::QuantumObject{DT2,StateOpType},
    tlist,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
    DT1,
    DT2,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    tlist = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    L_evo = _mesolve_make_L_QobjEvo(H, c_ops)
    check_dims(L_evo, ψ0)

    ρ0 = sparse_to_dense(_CType(ψ0), mat2vec(ket2dm(ψ0).data)) # Convert it to dense vector with complex element type
    L = L_evo.data

    progr = ProgressBar(length(tlist), enable = getVal(progress_bar))

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(tlist))
        e_ops_data = ()
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(tlist))
        e_ops_data = [_generate_mesolve_e_op(op) for op in e_ops]
        is_empty_e_ops = isempty(e_ops)
    end

    p = (
        e_ops = e_ops_data,
        expvals = expvals,
        progr = progr,
        times = tlist,
        Hdims = L_evo.dims,
        is_empty_e_ops = is_empty_e_ops,
        params...,
    )

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_mesolve_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2)

    tspan = (tlist[1], tlist[end])
    return ODEProblem{true,FullSpecialize}(L, ρ0, tspan, p; kwargs3...)
end

@doc raw"""
    mesolve(
        H::Union{AbstractQuantumObject{DT1,HOpType},Tuple},
        ψ0::QuantumObject{DT2,StateOpType},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::NamedTuple = NamedTuple(),
        progress_bar::Union{Val,Bool} = Val(true),
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
- `params`: `NamedTuple` of parameters to pass to the solver.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
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
    H::Union{AbstractQuantumObject{DT1,HOpType},Tuple},
    ψ0::QuantumObject{DT2,StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {
    DT1,
    DT2,
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
