export Propagator
export propagator

function propagatorProblem(
    H::AbstractQuantumObject{Operator},
    tspan::AbstractVector;
    params = NullParameters(),
    progress_bar::Union{Val, Bool} = Val(true),
    inplace::Union{Val, Bool} = Val(true),
    kwargs...
)
    H_evo = QobjEvo(H)

    p0 = qeye(H).data
    U = -1im*H_evo.data

    kwargs2 = _merge_saveat(tspan, nothing, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_se_me_kwargs(nothing, makeVal(progress_bar), tspan, kwargs2, SaveFuncSESolve)
    
    return ODEProblem{getVal(inplace),FullSpecialize}(U, p0, tspan, params; kwargs3...)
end

function propagatorProblem(
    H::AbstractQuantumObject{SuperOperator},
    tspan::AbstractVector;
    params = NullParameters(),
    progress_bar::Union{Val, Bool} = Val(true),
    inplace::Union{Val, Bool} = Val(true),
    kwargs...
)
    H_evo = QobjEvo(H)

    p0 = qeye(H).data
    U = H_evo.data

    kwargs2 = _merge_saveat(tspan, nothing, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_se_me_kwargs(nothing, makeVal(progress_bar), tspan, kwargs2, SaveFuncMESolve)
    
    return ODEProblem{getVal(inplace),FullSpecialize}(U, p0, tspan, params; kwargs3...)
end

struct Propagator{T<:Union{AbstractQuantumObject{Operator}, AbstractQuantumObject{SuperOperator}}}
    H::T
    times::AbstractArray{Float64}
    props::AbstractArray{A} where A <: AbstractQuantumObject
    tol::Float64
    solver_kwargs::Base.Pairs
    dims::Union{AbstractArray, Tuple}
    solver::OrdinaryDiffEqAlgorithm
    type 
end

function Base.show(io::IO, p::Propagator)
    println("Propagator: ")
    println("   dims: $(p.dims)")
    println("   issuper: $(issuper(p.H))")
    println("   size: $(size(p.H))")
    println("   Number of Saved Times: $(length(p.times))")
end

function propagator(
    H::Union{QobjEvo, QuantumObject};
    tol = 1e-6,
    solver = Tsit5(),
    kwargs...
    )
    return Propagator(H, [0.0], [qeye_like(H)], tol, kwargs, H.dims, solver, Operator())
end

function propagator(
    H::Union{QobjEvo, QuantumObject},
    c_ops :: Union{AbstractArray, Tuple};
    tol = 1e-6,
    solver = Tsit5(),
    kwargs...
    )
    L = liouvillian(H, c_ops)
    return Propagator(L, [0.0], [L], tol, kwargs, L.dims, solver, SuperOperator())
end

function _lookup_or_compute(p::Propagator, t::Real)
    """
    Get U(t) from cache or compute it.
    """
    idx = searchsorted(p.times, t)

    # Check if we found an exact match or close enough within tolerance
    if (idx.start <= length(p.times)) && (abs(t - p.times[idx.start]) <= p.tol)
        U = p.props[idx.start]
    elseif (idx.start > 1) && (abs(t - p.times[idx.start - 1]) <= p.tol)
        U = p.props[idx.start - 1]
    else
        t0 = (idx.start > 1) ? p.times[idx.start - 1] : 0.0
        U = _compute_propagator(p, t)
        insert!(p.props, idx.start, U)
        insert!(p.times, idx.start, t)
    end

    return U
end

function _compute_propagator(p::Propagator,
    t::Real;
    t0 = 0.0
    )
    prob = propagatorProblem(p.H, [t0, t]; p.solver_kwargs...)
    res = solve(prob, p.solver)
    U = QuantumObject(res.u[end], dims = p.dims, type = p.type)
    return U
end


"""
    (p::Propagator)(t; t0=0.0, just_interval=false)

Evaluate the propagator at time t.

- t: target time.
- t0: initial time. If nonzero:
  - when just_interval == true, integrate directly on [t0, t] to obtain U(t, t0).
  - otherwise, compute U(t) * inv(U(t0)) using cached values when available.
- just_interval: force a fresh integration on [t0, t] instead of composing from cache. This is useful if both t0 and t are large but t-t0 is small. 

Returns a QuantumObject with the same dims/type as p.H. Results are cached at requested
times (within tolerance p.tol) to avoid recomputation.
"""
function (p::Propagator)(t::Real; t0 = 0.0, just_interval = false)
    if just_interval
        return _compute_propagator(p, t; t0 = t0)
    end
    # We start by computing U_t0 if needed. This means that the computation for U will, at worst, start at t0 instead of 0.0.
    if t0 != 0.0
        U_t0 = (_lookup_or_compute(p, t0))
    end
    U = _lookup_or_compute(p, t)

    # We only do the multiplication if needed
    if t0 != 0.0
        U = U*inv(U_t0)
    end
    return U
end

function Base.size(p::Propagator)
    return size(p.H)
end

