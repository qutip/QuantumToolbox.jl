export mesolveProblem, mesolve, mesolve_map

_mesolve_make_L_QobjEvo(H::Union{QuantumObject,Nothing}, c_ops) = QobjEvo(liouvillian(H, c_ops); type = SuperOperator())
_mesolve_make_L_QobjEvo(H::Union{QuantumObjectEvolution,Tuple}, c_ops) = liouvillian(QobjEvo(H), c_ops)
_mesolve_make_L_QobjEvo(H::Nothing, c_ops::Nothing) = throw(ArgumentError("Both H and
c_ops are Nothing. You are probably running the wrong function."))

function _gen_mesolve_solution(sol, times, dimensions, isoperket::Val)
    if getVal(isoperket)
        ρt = map(ϕ -> QuantumObject(ϕ, type = OperatorKet(), dims = dimensions), sol.u)
    else
        ρt = map(ϕ -> QuantumObject(vec2mat(ϕ), type = Operator(), dims = dimensions), sol.u)
    end

    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionSol(
        times,
        sol.t,
        ρt,
        _get_expvals(sol, SaveFuncMESolve),
        sol.retcode,
        sol.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end

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

    return _gen_mesolve_solution(sol, prob.times, prob.dimensions, prob.kwargs.isoperket)
end

@doc raw"""
    mesolve_map(
        H::Union{AbstractQuantumObject,Tuple},
        ψ0::Union{QuantumObject,AbstractVector{<:QuantumObject}},
        tlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
        alg::OrdinaryDiffEqAlgorithm = Tsit5(),
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        params::Union{NullParameters,Tuple} = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Solve the master equation for multiple initial states and parameter sets using ensemble simulation.

This function computes the time evolution for all combinations (Cartesian product) of initial states and parameter sets, solving the Lindblad master equation (see [`mesolve`](@ref)):

```math
\frac{\partial \hat{\rho}(t)}{\partial t} = -i[\hat{H}, \hat{\rho}(t)] + \sum_n \mathcal{D}(\hat{C}_n) [\hat{\rho}(t)]
```

where

```math
\mathcal{D}(\hat{C}_n) [\hat{\rho}(t)] = \hat{C}_n \hat{\rho}(t) \hat{C}_n^\dagger - \frac{1}{2} \hat{C}_n^\dagger \hat{C}_n \hat{\rho}(t) - \frac{1}{2} \hat{\rho}(t) \hat{C}_n^\dagger \hat{C}_n
```

for each combination in the ensemble.

# Arguments

- `H`: Hamiltonian of the system ``\hat{H}``. It can be either a [`QuantumObject`](@ref), a [`QuantumObjectEvolution`](@ref), or a `Tuple` of operator-function pairs.
- `ψ0`: Initial state(s) of the system. Can be a single [`QuantumObject`](@ref) or a `Vector` of initial states. It can be either a [`Ket`](@ref), [`Operator`](@ref) or [`OperatorKet`](@ref).
- `tlist`: List of time points at which to save either the state or the expectation values of the system.
- `c_ops`: List of collapse operators ``\{\hat{C}_n\}_n``. It can be either a `Vector` or a `Tuple`.
- `alg`: The algorithm for the ODE solver. The default is `Tsit5()`.
- `ensemblealg`: Ensemble algorithm to use for parallel computation. Default is `EnsembleThreads()`.
- `e_ops`: List of operators for which to calculate expectation values. It can be either a `Vector` or a `Tuple`.
- `params`: A `Tuple` of parameter sets. Each element should be an `AbstractVector` representing the sweep range for that parameter. The function will solve for all combinations of initial states and parameter sets.
- `progress_bar`: Whether to show the progress bar. Using non-`Val` types might lead to type instabilities.
- `kwargs`: The keyword arguments for the ODEProblem.

# Notes

- The function returns an array of solutions with dimensions matching the Cartesian product of initial states and parameter sets.
- If `ψ0` is a vector of `m` states and `params = (p1, p2, ...)` where `p1` has length `n1`, `p2` has length `n2`, etc., the output will be of size `(m, n1, n2, ...)`.
- If `H` is an [`Operator`](@ref), `ψ0` is a [`Ket`](@ref) and `c_ops` is `Nothing`, the function will call [`sesolve_map`](@ref) instead.
- See [`mesolve`](@ref) for more details.

# Returns

- An array of [`TimeEvolutionSol`](@ref) objects with dimensions `(length(ψ0), length(params[1]), length(params[2]), ...)`.
"""
function mesolve_map(
    H::Union{AbstractQuantumObject{HOpType},Tuple},
    ψ0::AbstractVector{<:QuantumObject{StateOpType}},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params::Union{NullParameters,Tuple} = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {HOpType<:Union{Operator,SuperOperator},StateOpType<:Union{Ket,Operator,OperatorKet}}
    (isoper(H) && all(isket, ψ0) && isnothing(c_ops)) && return sesolve_map(
        H,
        ψ0,
        tlist;
        alg = alg,
        ensemblealg = ensemblealg,
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        kwargs...,
    )

    # mapping initial states and parameters
    # Convert to appropriate format based on state type
    ψ0_iter = map(ψ0) do state
        T = _complex_float_type(eltype(state))
        if isoperket(state)
            to_dense(T, copy(state.data))
        else
            to_dense(T, mat2vec(ket2dm(state).data))
        end
    end
    if params isa NullParameters
        iter = collect(Iterators.product(ψ0_iter, [params])) |> vec # convert nx1 Matrix into Vector
    else
        iter = collect(Iterators.product(ψ0_iter, params...))
    end

    # we disable the progress bar of the mesolveProblem because we use a global progress bar for all the trajectories
    prob = mesolveProblem(
        H,
        first(ψ0),
        tlist,
        c_ops;
        e_ops = e_ops,
        params = first(iter)[2:end],
        progress_bar = Val(false),
        kwargs...,
    )

    return mesolve_map(prob, iter, alg, ensemblealg; progress_bar = progress_bar)
end
mesolve_map(
    H::Union{AbstractQuantumObject{HOpType},Tuple},
    ψ0::QuantumObject{StateOpType},
    tlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    kwargs...,
) where {HOpType<:Union{Operator,SuperOperator},StateOpType<:Union{Ket,Operator,OperatorKet}} =
    mesolve_map(H, [ψ0], tlist, c_ops; kwargs...)

# this method is for advanced usage
# User can define their own iterator structure, prob_func and output_func
#   - `prob_func`: Function to use for generating the ODEProblem.
#   - `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `Progress` object, and the (optional) `RemoteChannel` object.
#
# Return: An array of TimeEvolutionSol objects with the size same as the given iter.
function mesolve_map(
    prob::TimeEvolutionProblem{<:ODEProblem},
    iter::AbstractArray,
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    ensemblealg::EnsembleAlgorithm = EnsembleThreads();
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
)
    # generate ensemble problem
    ntraj = length(iter)
    _prob_func = isnothing(prob_func) ? (prob, i, repeat) -> _se_me_map_prob_func(prob, i, repeat, iter) : prob_func
    _output_func =
        isnothing(output_func) ?
        _ensemble_dispatch_output_func(
            ensemblealg,
            progress_bar,
            ntraj,
            _standard_output_func;
            progr_desc = "[mesolve_map] ",
        ) : output_func
    ens_prob = TimeEvolutionProblem(
        EnsembleProblem(prob.prob, prob_func = _prob_func, output_func = _output_func[1], safetycopy = false),
        prob.times,
        prob.dimensions,
        (progr = _output_func[2], channel = _output_func[3], isoperket = prob.kwargs.isoperket),
    )

    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemblealg, ntraj)

    # handle solution and make it become an Array of TimeEvolutionSol
    sol_vec =
        [_gen_mesolve_solution(sol[:, i], prob.times, prob.dimensions, prob.kwargs.isoperket) for i in eachindex(sol)] # map is type unstable
    return reshape(sol_vec, size(iter))
end
