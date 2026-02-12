export Propagator, propagator
@doc raw"""
    Propagator{HT, PT, DT, KWT}

A callable struct representing a time-evolution propagator for a quantum system.

It lazily computes and caches propagators over requested time intervals. For time-independent Hamiltonians
([`QuantumObject`](@ref)), the propagator is computed via matrix exponentiation:

```math
\hat{U}(t, t_0) = e^{-i \hat{H} (t - t_0)}
```

For time-dependent Hamiltonians ([`QuantumObjectEvolution`](@ref)), the propagator is obtained by solving the
Schrödinger equation ([`sesolve`](@ref)) for [`Operator`](@ref) types, or the master equation ([`mesolve`](@ref))
for [`SuperOperator`](@ref) types, using an identity matrix as the initial condition.

Previously computed propagators for sub-intervals are reused automatically to avoid redundant computation.

# Fields

- `H`: The Hamiltonian or Liouvillian of the system.
- `props`: A dictionary mapping time intervals `[t0, t]` to their computed propagators.
- `dims`: The dimensions of the Hilbert space.
- `solver_kwargs`: Keyword arguments forwarded to the underlying solver ([`sesolve`](@ref) or [`mesolve`](@ref)).
- `max_saved`: Maximum number of propagators to cache.
- `threshold`: Numerical tolerance for matching stored time intervals.
- `remember_by_default`: Whether to cache newly computed propagators by default.
- `period`: If the system is periodic in time, the propagator will calculate intervals modulo this period.

# Usage

A `Propagator` is callable. Use `U(t; t0=0.0)` to obtain the propagator from `t0` to `t`, or `U([t0, t])` for
an interval. See [`propagator`](@ref) for construction.
"""
struct Propagator{
        HT <: Union{Operator, SuperOperator},
        PT <: AbstractQuantumObject,
        DT <: AbstractDimensions,
        KWT,
    }
    H::AbstractQuantumObject{HT}
    props::Dict{Vector, PT}
    dims::AbstractArray
    dimensions::DT
    solver_kwargs::KWT
    max_saved::Union{Integer, Float64}
    threshold::Float64
    remember_by_default::Bool
    period::Real
    isconstant::Bool
end


@doc raw"""
    propagator(
        H::AbstractQuantumObject{HOpType},
        t::Union{Nothing, Real} = nothing;
        t0 = 0.0,
        threshold::Float64 = 1e-9,
        max_saved::Union{Integer, Float64} = typemax(Int),
        remember_by_default::Bool = true,
        params = NullParameters(),
        progress_bar::Union{Val, Bool} = Val(true),
        inplace::Union{Val, Bool} = Val(true),
        kwargs...,
    )

Construct a [`Propagator`](@ref) object for the quantum system described by Hamiltonian (or Liouvillian) `H`.

If `t` is provided, the propagator from `t0` to `t` is immediately computed and cached. Otherwise, an empty
`Propagator` is returned, ready to be evaluated lazily at arbitrary times.

# Arguments

- `H`: Hamiltonian or Liouvillian of the system ``\hat{H}``. It can be a [`QuantumObject`](@ref) or a
  [`QuantumObjectEvolution`](@ref), with type [`Operator`](@ref) or [`SuperOperator`](@ref).
- `t`: Optional final time. If given, the propagator for the interval `[t0, t]` is computed immediately.
- `t0`: Initial time. Default is `0.0`.
- `threshold`: Numerical tolerance for matching cached time intervals. Default is `1e-9`.
- `period`: If the system is periodic in time, specify the period to automatically calculate intervals modulo this period. Default is `Inf` (no periodicity).
- `max_saved`: Maximum number of propagators to store in the cache. Can be an `Integer` or `Inf`. Default is `typemax(Int)`.
- `remember_by_default`: Whether to automatically cache computed propagators. Default is `true`.
- `params`: Parameters to pass to the underlying solver.
- `progress_bar`: Whether to show a progress bar during time evolution. Using non-`Val` types might lead to type instabilities.
- `inplace`: Whether to use inplace operations for the ODE solver. Default is `Val(true)`.
- `kwargs`: Additional keyword arguments forwarded to the solver ([`sesolve`](@ref) or [`mesolve`](@ref)).

# Returns

- `U::Propagator`: A callable propagator object. Use `U(t; t0=0.0)` or `U([t0, t])` to evaluate.

# Example

```julia
H = (ϵ / 2) * sigmaz() + (Ω / 2) * sigmax()
U = propagator(H)      # create lazy propagator
ψt = U(π) * basis(2, 0) # propagate |0⟩ to time π
```
"""
function propagator(
        H::AbstractQuantumObject{HOpType},
        t::Union{Nothing, Real} = nothing;
        t0 = 0.0,
        threshold::Float64 = 1.0e-9,
        period::Real = Inf,
        isconstant::Bool = false,
        max_saved::Union{Integer, Float64} = typemax(Int),
        remember_by_default::Bool = true,
        params = NullParameters(),
        progress_bar::Union{Val, Bool} = Val(true),
        inplace::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}}

    full_kwargs = (; params, progress_bar, inplace, kwargs...)

    if !(max_saved isa Integer) && max_saved != Inf
        max_saved = ceil(Int, max_saved)
        @warn "max_saved should be an Integer or Inf. Setting to $max_saved."
    end

    if !(H isa QobjEvo)
        isconstant = true
    end
    U = Propagator(H, Dict{Vector, AbstractQuantumObject}(), H.dims, H.dimensions, full_kwargs, max_saved, threshold, remember_by_default, period, isconstant)

    if t != nothing
        U(t; t0 = t0, remember = true)
    end
    return U
end


@doc raw"""
    (U::Propagator)(t; t0 = 0.0, remember = nothing, return_result = true, save_steps = true)

Evaluate the propagator `U` from time `t0` to time `t`.

This computes the time-evolution propagator ``\hat{U}(t, t_0)`` by combining any previously cached
sub-interval propagators with newly computed ones for uncovered gaps. The full interval `[t0, t]`
propagator is also cached separately when `remember` is enabled.

# Arguments

- `t`: The final time.
- `t0`: The initial time. Default is `0.0`.
- `remember`: Whether to cache newly computed propagators. If `nothing` (default), caching follows the
  `remember_by_default` setting of the [`Propagator`](@ref), subject to the `max_saved` limit.
- `return_result`: Whether to return the computed propagator. Default is `true`.
- `save_steps`: Whether to cache intermediate sub-interval propagators in addition to the full `[t0, t]`
  interval. Default is `true`. Set to `false` to only cache the composite result.

# Returns

- The propagator as an `AbstractQuantumObject` (if `return_result` is `true`).
"""
function (U::Propagator)(t; t0 = 0.0, remember::Union{Nothing, Bool} = nothing, return_result = true, save_steps = true)
    ΔT = abs(t - t0)
    if U.period != Inf
        t = mod(t, U.period)
        t0 = mod(t0, U.period)
    end
    intervals = _get_intervals_for_range(collect(keys(U.props)), [t0, t]; threshold = U.threshold)

    prop = qeye_like(U.H)
    if prop isa QobjEvo
        prop = prop(0.0)
    end

    all_intervals = vcat(intervals.usable, intervals.to_compute)
    sort!(all_intervals, by = first)

    for interval in all_intervals
        temp_prop = _propagator_compute_or_look_up(U, interval)
        if (remember === nothing ? (U.remember_by_default && length(U.props) < U.max_saved) : remember) && !(interval in keys(U.props)) && save_steps
            if length(U.props) >= U.max_saved
                @warn "Maximum number of stored propagators reached, save is being forced because 'remember' is set to true."
            end
            U.props[interval] = temp_prop
        end
        prop = temp_prop * prop
    end

    if U.period != Inf && (t - t0) > U.period
        prop = prop^(ΔT / U.period)
    end

    if (remember === nothing ? (U.remember_by_default && length(U.props) < U.max_saved) : remember) && !([t0, t] in keys(U.props))
        if length(U.props) >= U.max_saved
            @warn "Maximum number of stored propagators reached, save is being forced because 'remember' is set to true."
        end
        U.props[[t0, t]] = prop
    end

    if return_result
        return prop
    end
end


function (U::Propagator)(interval::Vector; kwargs...)
    return U(interval[2]; t0 = interval[1], kwargs...)
end


"""
    _propagator_compute_or_look_up(U::Propagator{HT}, interval) where HT

Look up a cached propagator for the given `interval`, or compute it if not found.

For time-independent Hamiltonians ([`QuantumObject`](@ref)), the interval is shifted to `[0, Δt]` since the
propagator depends only on the duration. For time-dependent Hamiltonians ([`QuantumObjectEvolution`](@ref)),
the propagator is computed via [`sesolve`](@ref) (for [`Operator`](@ref)) or [`mesolve`](@ref) (for
[`SuperOperator`](@ref)) using an identity matrix as the initial state.
"""
function _propagator_compute_or_look_up(U::Propagator{HT}, interval) where {HT <: Union{Operator, SuperOperator}}
    if U.isconstant
        interval = [0.0, interval[2] - interval[1]]
    end

    if interval in keys(U.props)
        return U.props[interval]
    else
        if U.isconstant
            if HT <: Operator
                return exp(-1im * U.H * (interval[2] - interval[1]))
            else
                return exp(U.H * (interval[2] - interval[1]))
            end
        end

        if HT <: Operator
            return sesolve(U.H, qeye_like(U.H)(0.0)::QuantumObject{Operator}, interval; saveat = [interval[2]], U.solver_kwargs...).states[end]
        else
            return mesolve(U.H, qeye_like(U.H)(0.0)::QuantumObject{SuperOperator}, interval; saveat = [interval[2]], U.solver_kwargs...).states[end]
        end
    end
end


"""
    _get_intervals_for_range(stored_intervals, target_interval; threshold=1e-9)

Decompose `target_interval = [a, b]` into sub-intervals by reusing `stored_intervals` where possible.

Returns a `NamedTuple` with:
- `usable`: Stored intervals fully contained within `[a, b]` (within `threshold` tolerance).
- `to_compute`: Gap intervals not covered by any stored interval that still need to be computed.
"""
function _get_intervals_for_range(stored_intervals::AbstractVector{T}, target_interval::Vector; threshold = 1.0e-9) where {T <: Vector}
    a, b = target_interval

    # Find stored intervals that are fully contained within target range (with fuzzy boundaries)
    usable = [[s, e] for (s, e) in stored_intervals if s >= a - threshold && e <= b + threshold]
    sort!(usable, by = first)

    # Merge usable intervals to find coverage
    merged = Vector[]
    for (s, e) in usable
        if isempty(merged) || s > merged[end][2] + threshold
            push!(merged, [s, e])
        else
            merged[end] = (merged[end][1], max(merged[end][2], e))
        end
    end

    # Find gaps that need to be computed
    to_compute = Vector[]
    current = a
    for (s, e) in merged
        if current < s - threshold
            push!(to_compute, [current, s])
        end
        current = max(current, e)
    end
    if current < b - threshold
        push!(to_compute, [current, b])
    end
    return (usable = usable, to_compute = to_compute)
end


function Base.show(io::IO, U::Propagator)
    saved_times = String[]
    times = collect(keys(U.props))
    for i in 1:length(times)
        push!(saved_times, "\n  $(times[i][1]) -> $(times[i][2])")
    end
    if length(saved_times) == 0
        saved_times = ["None"]
    end
    return println(
        io,
        "\nPropagator: H Type=",
        U.H.type,
        "   dims=",
        _get_dims_string(U.dimensions),
        "   size=",
        size(U),
        "\nperiod: ",
        U.period,
        "\nSaved Propagators: ",
        saved_times...,
        "\nMemory Usage: ",
        Base.format_bytes(Base.summarysize(U.props))
    )
end

function Base.size(U::Propagator)
    return size(U.H)
end

function Base.length(U::Propagator)
    return length(U.H)
end

function Base.pop!(U::Propagator, interval::Vector)
    if U.period != Inf
        interval = mod.(interval, U.period)
    end
    if interval in keys(U.props)
        return pop!(U.props, interval)
    else
        @warn "Interval $interval not found in cache. No propagator removed."
        return nothing
    end
end