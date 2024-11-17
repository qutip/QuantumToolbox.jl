export sesolveProblem, sesolve

function _merge_sesolve_kwargs_with_callback(kwargs, cb)
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(kwargs.callback, cb),)) :
        merge(kwargs, (callback = cb,))

    return kwargs2
end

# Multiple dispatch depending on the progress_bar and e_ops types
function _generate_sesolve_kwargs(e_ops, progress_bar, tlist, kwargs)
    cb = _generate_sesolve_callback(e_ops, tlist, progress_bar)
    return _merge_sesolve_kwargs_with_callback(kwargs, cb)
end
_generate_sesolve_kwargs(e_ops::Nothing, progress_bar::Val{false}, tlist, kwargs) = kwargs

_sesolve_make_U_QobjEvo(H::QuantumObjectEvolution{<:MatrixOperator}) =
    QobjEvo(MatrixOperator(-1im * H.data.A), dims = H.dims, type = Operator)
_sesolve_make_U_QobjEvo(H) = QobjEvo(H, -1im)

@doc raw"""
    sesolveProblem(
        H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{DT2,KetQuantumObject},
        tlist::AbstractVector;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NamedTuple, AbstractVector, TimeEvolutionParameters} = eltype(ψ0)[],
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
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` or `AbstractVector` of parameters to pass to the solver. For more advanced usage, you can use the [`TimeEvolutionParameters`](@ref) struct.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use the inplace version of the ODEProblem. The default is `Val(true)`.
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
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NamedTuple,AbstractVector,TimeEvolutionParameters} = eltype(ψ0)[],
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in QuantumToolbox."))

    tlist = convert(Vector{_FType(ψ0)}, tlist) # Convert it to support GPUs and avoid type instabilities for OrdinaryDiffEq.jl

    H_evo = _sesolve_make_U_QobjEvo(H) # Multiply by -i
    isoper(H_evo) || throw(ArgumentError("The Hamiltonian must be an Operator."))
    check_dims(H_evo, ψ0)

    ψ0 = sparse_to_dense(_CType(ψ0), get_data(ψ0)) # Convert it to dense vector with complex element type
    U = H_evo.data

    if e_ops isa Nothing
        expvals = Array{ComplexF64}(undef, 0, length(tlist))
        is_empty_e_ops = true
    else
        expvals = Array{ComplexF64}(undef, length(e_ops), length(tlist))
        is_empty_e_ops = isempty(e_ops)
    end

    if params isa TimeEvolutionParameters
        (e_ops isa Nothing) || throw(
            ArgumentError(
                "The parameter `params` cannot be a TimeEvolutionParameters object when `e_ops` is not Nothing",
            ),
        )
    end

    p = params isa TimeEvolutionParameters ? params : TimeEvolutionParameters(params, expvals)

    saveat = is_empty_e_ops ? tlist : [tlist[end]]
    default_values = (DEFAULT_ODE_SOLVER_OPTIONS..., saveat = saveat)
    kwargs2 = merge(default_values, kwargs)
    kwargs3 = _generate_sesolve_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs2)

    tspan = (tlist[1], tlist[end])
    prob = ODEProblem{getVal(inplace),FullSpecialize}(U, ψ0, tspan, p; kwargs3...)

    return TimeEvolutionProblem(prob, tlist, H_evo.dims)
end

@doc raw"""
    sesolve(
        H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
        ψ0::QuantumObject{DT2,KetQuantumObject},
        tlist::AbstractVector;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NamedTuple, AbstractVector} = eltype(ψ0)[],
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
- `tlist`: List of times at which to save either the state or the expectation values of the system.
- `alg`: The algorithm for the ODE solver. The default is `Tsit5()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: `NamedTuple` or `AbstractVector` of parameters to pass to the solver.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use the inplace version of the ODEProblem. The default is `Val(true)`.
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
    H::Union{AbstractQuantumObject{DT1,OperatorQuantumObject},Tuple},
    ψ0::QuantumObject{DT2,KetQuantumObject},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NamedTuple,AbstractVector} = eltype(ψ0)[],
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {DT1,DT2}
    prob = sesolveProblem(
        H,
        ψ0,
        tlist;
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )

    return sesolve(prob, alg)
end

function sesolve(prob::TimeEvolutionProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5())
    sol = solve(prob.prob, alg)

    ψt = map(ϕ -> QuantumObject(ϕ, type = Ket, dims = prob.dims), sol.u)

    return TimeEvolutionSol(
        prob.times,
        ψt,
        sol.prob.p.expvals,
        sol.retcode,
        sol.alg,
        sol.prob.kwargs[:abstol],
        sol.prob.kwargs[:reltol],
    )
end
