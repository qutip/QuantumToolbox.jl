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

mesolve_ti_dudt!(du, u, p, t) = mul!(du, p.L, u)
function mesolve_td_dudt!(du, u, p, t)
    mul!(du, p.L, u)
    L_t = p.H_t(t, p)
    return mul!(du, L_t, u, 1, 1)
end

"""
    mesolveProblem(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector=[];
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Bool=true,
        kwargs...)

Generates the ODEProblem for the master equation time evolution of an open quantum system.

# Arguments

  - `H::QuantumObject`: The Hamiltonian or the Liouvillian of the system.
  - `ψ0::QuantumObject`: The initial state of the system.
  - `t_l::AbstractVector`: The time list of the evolution.
  - `c_ops::AbstractVector=[]`: The list of the collapse operators.
  - `alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm=Tsit5()`: The algorithm used for the time evolution.
  - `e_ops::AbstractVector=[]`: The list of the operators for which the expectation values are calculated.
  - `H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing`: The time-dependent Hamiltonian or Liouvillian.
  - `params::NamedTuple=NamedTuple()`: The parameters of the time evolution.
  - `progress_bar::Bool=true`: Whether to show the progress bar.
  - `kwargs...`: The keyword arguments for the ODEProblem.

# Returns

  - `prob::ODEProblem`: The ODEProblem for the master equation time evolution.
"""
function mesolveProblem(
    H::QuantumObject{MT1,HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l,
    c_ops::Vector{QuantumObject{Tc,COpType}} = QuantumObject{MT1,HOpType}[];
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{Te,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Bool = true,
    kwargs...,
) where {
    MT1<:AbstractMatrix,
    T2,
    Tc<:AbstractMatrix,
    Te<:AbstractMatrix,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
    COpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
}
    H.dims != ψ0.dims && throw(ErrorException("The two operators don't have the same Hilbert dimension."))

    is_time_dependent = !(H_t === nothing)

    ρ0 = mat2vec(ket2dm(ψ0).data)
    L = liouvillian(H, c_ops).data

    progr = ProgressBar(length(t_l), enable = progress_bar)
    expvals = Array{ComplexF64}(undef, length(e_ops), length(t_l))
    e_ops2 = @. mat2vec(adjoint(get_data(e_ops)))

    p = (
        L = L,
        progr = progr,
        Hdims = H.dims,
        e_ops = e_ops2,
        expvals = expvals,
        H_t = H_t,
        is_empty_e_ops = isempty(e_ops),
        params...,
    )

    default_values = (abstol = 1e-7, reltol = 1e-5, saveat = [t_l[end]])
    kwargs2 = merge(default_values, kwargs)
    if !isempty(e_ops) || progress_bar
        cb1 = PresetTimeCallback(t_l, _save_func_mesolve, save_positions = (false, false))
        kwargs2 =
            haskey(kwargs, :callback) ? merge(kwargs2, (callback = CallbackSet(kwargs2.callback, cb1),)) :
            merge(kwargs2, (callback = cb1,))
    end

    tspan = (t_l[1], t_l[end])
    return _mesolveProblem(L, ρ0, tspan, alg, Val(is_time_dependent), p; kwargs2...)
end

function _mesolveProblem(
    L::AbstractMatrix{<:T1},
    ρ0::AbstractVector{<:T2},
    tspan::Tuple,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    is_time_dependent::Val{false},
    p;
    kwargs...,
) where {T1,T2}
    return ODEProblem{true,SciMLBase.FullSpecialize}(mesolve_ti_dudt!, ρ0, tspan, p; kwargs...)
end

function _mesolveProblem(
    L::AbstractMatrix{<:T1},
    ρ0::AbstractVector{<:T2},
    tspan::Tuple,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
    is_time_dependent::Val{true},
    p;
    kwargs...,
) where {T1,T2}
    return ODEProblem{true,SciMLBase.FullSpecialize}(mesolve_td_dudt!, ρ0, tspan, p; kwargs...)
end

"""
    mesolve(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector, c_ops::AbstractVector=[];
        alg::OrdinaryDiffEqAlgorithm=Tsit5(),
        e_ops::AbstractVector=[],
        H_t::Union{Nothing,Function,TimeDependentOperatorSum}=nothing,
        params::NamedTuple=NamedTuple(),
        progress_bar::Bool=true,
        kwargs...)

Time evolution of an open quantum system using master equation.

# Arguments

  - `H::QuantumObject`: Hamiltonian of Liouvillian of the system.
  - `ψ0::QuantumObject`: Initial state of the system.
  - `t_l::AbstractVector`: List of times at which to save the state of the system.
  - `c_ops::AbstractVector`: List of collapse operators.
  - `alg::OrdinaryDiffEqAlgorithm`: Algorithm to use for the time evolution.
  - `e_ops::AbstractVector`: List of operators for which to calculate expectation values.
  - `H_t::Union{Nothing,Function,TimeDependentOperatorSum}`: Time-dependent part of the Hamiltonian.
  - `params::NamedTuple`: Named Tuple of parameters to pass to the solver.
  - `progress_bar::Bool`: Whether to show the progress bar.
  - `kwargs...`: Additional keyword arguments to pass to the solver.

# Returns

  - `sol::TimeEvolutionSol`: The solution of the time evolution.
"""
function mesolve(
    H::QuantumObject{MT1,HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector,
    c_ops::Vector{QuantumObject{Tc,COpType}} = QuantumObject{MT1,HOpType}[];
    alg::OrdinaryDiffEqAlgorithm = Tsit5(),
    e_ops::Vector{QuantumObject{Te,OperatorQuantumObject}} = QuantumObject{MT1,OperatorQuantumObject}[],
    H_t::Union{Nothing,Function,TimeDependentOperatorSum} = nothing,
    params::NamedTuple = NamedTuple(),
    progress_bar::Bool = true,
    kwargs...,
) where {
    MT1<:AbstractMatrix,
    T2,
    Tc<:AbstractMatrix,
    Te<:AbstractMatrix,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
    COpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
}
    prob = mesolveProblem(
        H,
        ψ0,
        t_l,
        c_ops;
        alg = alg,
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        progress_bar = progress_bar,
        kwargs...,
    )

    return mesolve(prob, alg; kwargs...)
end

function mesolve(prob::ODEProblem, alg::OrdinaryDiffEqAlgorithm = Tsit5(); kwargs...)
    sol = solve(prob, alg)

    return _mesolve_sol(sol; kwargs...)
end

function _mesolve_sol(sol; kwargs...)
    Hdims = sol.prob.p.Hdims
    ρt = !haskey(kwargs, :save_idxs) ? map(ϕ -> QuantumObject(vec2mat(ϕ), dims = Hdims), sol.u) : sol.u

    return TimeEvolutionSol(sol.t, ρt, sol.prob.p.expvals)
end
