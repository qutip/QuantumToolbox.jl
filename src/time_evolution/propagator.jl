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
    H_evo = _sesolve_make_U_QobjEvo(H)
    U = H_evo.data

    p0 = one(H_evo(0)).data

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
    L_evo = _mesolve_make_L_QobjEvo(H, [])
    U = L_evo.data

    p0 = one(L_evo(0)).data

    kwargs2 = _merge_saveat(tspan, nothing, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _generate_se_me_kwargs(nothing, makeVal(progress_bar), tspan, kwargs2, SaveFuncMESolve)
    
    return ODEProblem{getVal(inplace),FullSpecialize}(U, p0, tspan, params; kwargs3...)
end

@doc raw"""
Propagator{T}

Container for time-evolution propagators built from `H` (a Hamiltonian or a Louivillian superoperator).
This immutable struct groups the generator, time list and solver settings together with any cached
propagators so that previous results can be reused where applicable.

# Parameters
- `H::T`
    The generator of dynamics. `T` is typically an `AbstractQuantumObject{Operator}` for closed,
    unitary dynamics (Hamiltonian) or an `AbstractQuantumObject{SuperOperator}` for open-system
    dynamics (Liouvillian or other superoperator).
- `times::AbstractArray{Float64}`
    A monotonic array of time points at which propagators are requested or saved. Values are
    interpreted in the same time units used to build `H`.
- `props::AbstractArray{A} where A <: AbstractQuantumObject`
    An array holding computed propagators corresponding to `times`. Elements are quantum objects
    (operators or superoperators) matching the expected output of propagating with `H`. Both `times` and `props` are initialized with the identity operator/superoperator at t=0.
- `tol::Float64`
    This is the time tolerance used to determine if a previous result should be reused. If \$t - t_{\text{cached}} < \text{tol}\$, where $t_\text{cached}\$ is the closed cached
    time to \$t\$, then the propagator at \$t_\text{cached}\$ will be reused.
- `solver_kwargs::Base.Pairs`
    Additional keyword arguments forwarded to the ODE integrator (for example `abstol`, `reltol`,
    `saveat`, or solver-specific options). Stored as a `Base.Pairs` collection for convenient
    dispatch into solver APIs.
- `dims::Union{AbstractArray, Tuple}`
    Dimension metadata of `H`.
- `solver::OrdinaryDiffEqAlgorithm`
    The chosen ODE algorithm (a concrete algorithm type from OrdinaryDiffEq, the default is `Tsit5`) that will be used to integrate the generator to obtain
    time-ordered evolution.
- `type`
    Whether the propagator is an operator or superoperator. 

# Initialization
A propagator is initialized by calling the `propagator` function.

# Calling
If `U` is a propagator object, to get the propagator between t0 and t, call `U(t, t0=t0; just_interval::Bool=false)`.
    - t0 is optional, if left out then t0=0.0. 
    - `just_interval` determines how the propagator from t0 to t is calculated. 
        - if `just_interval = false`, then the propagator is calculated by getting U1(t0,0.0) and U2(t, 0.0) and setting
        ```math
        U = U2U1^{-1}
        ```
        This is useful if U1 is alread cached or one plans to reuse U1 or U2. 
        - if `just_interval=true` then the propagator is calulated by directly integrating from t0 to t. This approach does't save 
        anything but is useful if t-t0 >> t0 as it avoids the long integration required to get t0. 
"""
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

@doc raw"""
propagator(H::Union{QobjEvo, QuantumObject}; tol=1e-6, solver=Tsit5(), kwargs...)

Construct and return a Propagator object prepared for time evolution using the
given Hamiltonian or time-dependent operator.

# Arguments
- H::Union{QobjEvo, QuantumObject}
    The Hamiltonian or generator for the dynamics. Can be a time-independent
    QuantumObject or a time-dependent QobjEvo describing H(t).
- c_ops::Union{AbstractArray, Tuple} (Optional)
    - if no c_ops are given, then this is treated as a lossless propagator and the propagator itself will be of type `AbstractQuantumObject{Operator}`. If c_ops are given, 
    then the propagator will be of type `AbstractQuantumObject{SuperOperator}` and the save `H` will be the Louivillian.
- tol::Real=1e-6
    Absolute/relative tolerance used by the ODE integrator. This value is
    forwarded to the integrator/Propagator to control numerical accuracy.
- solver
    The ODE solver algorithm to use (default: `Tsit5()`); any solver object
    compatible with the underlying ODE interface is accepted.
- kwargs...
    Additional keyword arguments are forwarded to the Propagator constructor
    (e.g., options for step control, callback functions, or integration
    settings supported by the backend).

# Notes
- The function does not perform propagation itself; it only constructs and
  configures the Propagator object which can then be used to evolve states or
  operators.
- For time-dependent `H` provide a `QobjEvo` with an appropriate time
  dependence. For time-independent dynamics provide a `QuantumObject`.

# Returns
- `Propagator`
    An initialized Propagator instance set up with
    - initial time [0.0],
    - initial operator equal to the identity on the Hilbert space of `H`
      (created with `one(H)`),
    - the requested tolerance and solver,
    - dimensions inferred from `H`, and
    - any extra options passed via `kwargs`.
"""
function propagator(
    H::Union{QobjEvo, QuantumObject};
    tol = 1e-6,
    solver = Tsit5(),
    kwargs...
    )
    H_evo = QobjEvo(H)
    return Propagator(H, [0.0], [one(H_evo(0))], tol, kwargs, H.dims, solver, Operator())
end

function propagator(
    H::Union{QobjEvo, QuantumObject},
    c_ops :: Union{AbstractArray, Tuple};
    tol = 1e-6,
    solver = Tsit5(),
    kwargs...
    )
    L = liouvillian(H, c_ops)
    L_evo = QobjEvo(L)
    return Propagator(L, [0.0], [one(L_evo(0))], tol, kwargs, L.dims, solver, SuperOperator())
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


function (p::Propagator)(t::Real, t0 = 0.0; just_interval = false)
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

