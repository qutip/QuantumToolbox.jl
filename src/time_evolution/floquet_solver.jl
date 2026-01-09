using Revise
# script helper functions
function _to_period_interval(tlist::AbstractVector, T::Real)
    # function maps all elements ``t`` in `tlist` outside the interval ``[0, T)`` to an equivalent
    # time ``\tau`` such that ``mod(t, T) = \tau``
    if !isempty(tlist)
        tlist = mod.(tlist, T)
        unique!(tlist)
        sort!(tlist)
    end
    return tlist
end



@doc raw"""
   struct FloquetBasis

Julia struct containing propagators, quasienergies, and Floquet states for a system with a ``T``-periodic Hamiltonain.

# Fields:
- `H::AbstractQuantumObject`: T-periodic Hamiltonian.
- `T<:Real`: Hamiltonian period such that ``\hat{H}(T) = \hat{H}(0)``
- `tlist::TE`: Time array fed to `sesolve` to compute propagators and Floquet states. First and final elements are always `0, T`. All elements lie in range `[0,T]`, see notes for behavior when field is set to an array with points outside this range.
- `precompute::TE`: Times for which the micromotion propagator and Floquet modes are precomputed. All elements  lie in range `[0,T]`. See notes for behavior when field is set to an array with points outside this range.
- `U_T::Qobj`: System propagator at time ``T``
- `Ulist::AbstractVector{Qobj}`: List of system propagators at times `tlist`
- `equasi::TE`: Time-independent quasienergies
"""
struct FloquetBasis{
    TT<:AbstractVector{<:Real},
    TE<:AbstractVector{<:Real},
    TQ<:AbstractVector{<:AbstractQuantumObject},
}
    H::AbstractQuantumObject
    T::Real
    tlist::TT
    precompute::TT
    U_T::Qobj
    Ulist::TQ
    equasi::TE

    @doc raw"""
        FloquetBasis(H::AbstractQuantumObject, T::Real, tlist::AbstractVector{Real}, precompute::Bool = true; kwargs::Dict = Dict())

        DOCSTRING

        # Arguments:
        - `H`: Time-dependent system Hamiltonian.
        - `T`: Hamiltonian period such that ``\hat{H}(T+\tau_0) = \hat{H}(\tau_0)``.
        - `tlist`: Time vector to use internally in sesolve to calculate the period-propagator and quasienergies.
        - `precompute`: Time vector containing points ``t`` at which to store the system propagator ``U(t)``.
        - `kwargs`: Additional keyword arguments to pass to ssesolve.

        # Notes:
        - If `tlist` or `precompute` contain elements outside the interval ``[0,T]``, a new time vector will be produced with all times ``t_k`` not in the interval mapped to an equivalent time ``\tau_k`` in the interval such that ``\hat{H}(t_k) = \hat{H}(\tau_k)``.
        - If the first and final elements of `tlist` are not 0 and T, then 0 and T will be prepended and appended to `tlist`.
        - If all elements of `precompute` are not in `tlist`, `tlist` will be set to the union of the two vectors.

        # Returns:
        - `fbasis::FloquetBasis`: FloquetBasis object for the system evolving under the ``T``-periodic Hamiltonian `H`.
        """
    function FloquetBasis(
        H::AbstractQuantumObject,
        T::Real,
        tlist::TT,
        precompute::TT;
        kwargs::Dict=Dict()
        ) where {TT<:AbstractVector{<:Real}}
        if T<=0
            throw(
                ArgumentError(
                    "`T` must be a nonzero positive real number"
                              ))
        else
            tlist, precompute = _to_period_interval.([tlist, precompute], T) # enforce that all timepoints lie in interval [0,T)
            tlist = union(tlist, precompute) # ensure all times in precompute are in tlist
            tlist, precompute = [unique((0, tlist..., T)), unique((precompute..., T))] # ensure that period-propagator is calculated
        end
        kwargs[:saveat] = precompute
        Ulist = sesolve(H, qeye_like(H), tlist, kwargs=kwargs).states
        U_T = pop!(Ulist)
        period_phases = eigenenergies(U_T)
        equasi = angle.(period_phases) ./ T
        new{typeof(tlist), typeof(equasi), typeof(Ulist)}(H, T, tlist, precompute,  U_T, Ulist, equasi)
    end
end

@doc raw"""
    FloquetBasis(H::AbstractQuantumObject, T::Real; kwargs::Dict = Dict())

DOCSTRING

# Arguments:
- `H`: Time-dependent system Hamiltonian.
- `T`: Hamiltonian period such that ``\hat{H}(T+\tau_0) = \hat{H}(\tau_0)``.
- `precompute`: If true, resulting FloquetBasis object will store precomputed propagators for all times in `range(start:0, stop:T, length:101)`. If false, only the final period propagator ``U(T)`` will be stored. Default is `true`.
- `kwargs`: Additional keyword arguments to pass to ssesolve.

# Notes
- Calling `FloquetBasis` without providing a time vector will create a FloquetBasis object with default `tlist=range(start:0, stop:T, length:101)`.

# Returns:
- `fbasis::FloquetBasis`: Floquet basis object for the system evolving under the time-dependent Hamiltonian `H`.
"""
function FloquetBasis(H::AbstractQuantumObject, T::Real, precompute::Bool=true; kwargs::Dict=Dict())
    tlist = range(0, T, 101)
    precompute = Float64[]
    return FloquetBasis(H, T, tlist, precompute, kwargs=kwargs)
end

@doc raw"""
    FloquetBasis(H::AbstractQuantumObject, T::Real, tlist::AbstractVector{Real}, precompute::Bool = true; kwargs::Dict = Dict())

DOCSTRING

# Arguments:
- `H`: Time-dependent system Hamiltonian.
- `T`: Hamiltonian period such that ``\hat{H}(T+\tau_0) = \hat{H}(\tau_0)``.
- `tlist`: Time vector to use internally in sesolve to calculate Period and intra-period propagators (if `precompute` is not `false`)
- `precompute`: If true, resulting FloquetBasis object will store precomputed propagators for all times in `tlist`. If false, only the final period propagator ``U(T)`` will be stored. Default is `true`.
- `kwargs`: Additional keyword arguments to pass to ssesolve.

# Notes:
- If `tlist` contains elements outside the interval ``[0,T]``, a new `tlist` will be produced with all times ``t_k`` not in the interval mapped to an equivalent time ``\tau_k`` in the interval such that ``\hat{H}(t_k) = \hat{H}(\tau_k)``.
- If the first and final elements of `tlist` are not 0 and T, then 0 and T will be prepended and appended to `tlist`.

# Returns:
- `fbasis::FloquetBasis`: Floquet basis object for the system evolving under the time-dependent Hamiltonian `H`.
"""
function FloquetBasis(
    H::AbstractQuantumObject,
    T::Real,
    tlist::AbstractVector{Real},
    precompute::Bool=true;
    kwargs::Dict=Dict())
    if precompute
        return FloquetBasis(H, T, tlist, tlist, kwargs=kwargs)
    else
        return FloquetBasis(H, T, tlist, Float64[], kwargs=kwargs)
    end
end
