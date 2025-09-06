export EnsembleTimeEvolutionProblem
"""
    EnsembleTimeEvolutionProblem{PT<:TimeEvolutionProblem, PF<:Function}

A structure representing an ensemble time evolution problem for quantum systems.

# Fields

  - `prob::PT`: The base time evolution problem.
  - `func::PF`: A function used to modify or sample parameters for each trajectory in the ensemble.
  - `iterate_params::Bool`: If `true`, parameters are iterated for each trajectory; otherwise, the same parameters are used.
  - `full_iterator::AbstractArray`: An array containing all parameter sets or states to be used in the ensemble.
  - `n_states::Int`: The number of initial states.
  - `trajectories::Int`: The total number of trajectories to simulate.

# Usage

This is used when setting up ensemble sesolve problems, useful for simulating multiple quantum states or parameter sets in parallel.

Example:

```julia
H = 2 * π * 0.1 * sigmax()
ψ0 = basis(2, 0) # spin-up
tlist = LinRange(0.0, 100.0, 100)

ψs = [ψ0, basis(2, 1)] # spin-up and spin-down

params = collect(Iterators.product([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]))
res = sesolve(H, ψs, tlist; params = params, iterate_params = true, alg = Tsit5(), progress_bar = false);
```
"""
struct EnsembleTimeEvolutionProblem{
    PT<:TimeEvolutionProblem,
    PF<:Function,
    X<:Vector{T} where T<:QuantumObject{Ket},
    Y<:AbstractArray,
}
    prob::PT
    func::PF
    states::X
    params::Y
    problem_dims::Tuple
    trajectories::Int
end
function EnsembleTimeEvolutionProblem(
    prob::PT,
    states::Vector{T},
    params::AbstractArray = [NullParameters()],
) where {PT<:TimeEvolutionProblem,T<:QuantumObject{Ket}}
    problem_dims = (length(states), length(params))

    function ensemble_func(prob, i, repeat)
        state_id = mod1(i, problem_dims[1])
        param_id = div(i - 1, problem_dims[1]) + 1
        return remake(prob, u0 = states[state_id].data, p = params[param_id])
    end
    trajectories = prod(problem_dims)
    return EnsembleTimeEvolutionProblem(prob, ensemble_func, states, params, problem_dims, trajectories)
end

function sesolve(
    prob::EnsembleTimeEvolutionProblem,
    alg::OrdinaryDiffEqAlgorithm = Tsit5();
    backend = EnsembleThreads(),
)
    ensemble_prob = EnsembleProblem(prob.prob.prob, prob_func = prob.func)
    sols = solve(ensemble_prob, alg, backend, trajectories = prob.trajectories)

    to_return = Array{TimeEvolutionSol}(undef, prob.problem_dims)
    for i in 1:length(sols)
        ψt = map(ϕ -> QuantumObject(ϕ, type = Ket(), dims = prob.prob.dimensions), sols[i].u)
        sol = TimeEvolutionSol(
            prob.prob.times,
            sols[i].t,
            ψt,
            _get_expvals(sols[i], SaveFuncSESolve),
            sols[i].retcode,
            sols[i].alg,
            sols[i].prob.kwargs[:abstol],
            sols[i].prob.kwargs[:reltol],
        )
        to_return[CartesianIndices(to_return)[i]] = sol
    end

    return to_return
end

function sesolve(
    H::Union{AbstractQuantumObject{Operator},Tuple},
    ψ0s::Vector{T},
    tlist::AbstractVector;
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(false),
    inplace::Union{Val,Bool} = Val(true),
    backend = EnsembleThreads(),
    kwargs...,
) where {T<:QuantumObject{Ket}}
    prob_init = sesolveProblem(
        H,
        ψ0s[1],
        tlist;
        e_ops = e_ops,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )

    trajectories = length(ψ0s)

    # function ensemble_func(prob, i, repeat)
    #     return remake(prob, u0 = ψ0s[i].data)
    # end

    # ensemble_prob = EnsembleTimeEvolutionProblem(prob_init, ensemble_func, ψ0s, trajectories)
    ensemble_prob = EnsembleTimeEvolutionProblem(prob_init, ψ0s, [params])
    return sesolve(ensemble_prob, alg; backend = backend)
end
