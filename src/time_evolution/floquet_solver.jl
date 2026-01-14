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


struct FloquetEvolutionSol{
    TT1<:AbstractVector{<:Real},
    TT2<:AbstractVector{<:Real},
    TS<:AbstractVector{<:AbstractQuantumObject},
    TE<:Union{AbstractMatrix{<:Real}, Nothing},
    AlgT<:AbstractODEAlgorithm,
    TolT<:Real
} <: TimeEvolutionSol
    times::TT1
    times_states::TT2
    states::TS
    expect::TE
    alg::AlgT
    abstol::TolT
    reltol::TolT
end



@doc raw"""
   struct FloquetBasis

Julia struct containing propagators, quasienergies, and Floquet states for a system with a ``T``-periodic Hamiltonain.

# Fields:
- `H::AbstractQuantumObject`: T-periodic Hamiltonian.
- `T<:Real`: Hamiltonian period such that ``\hat{H}(T) = \hat{H}(0)``
- `tlist::TT`: Time array fed to `sesolve` to compute propagators and Floquet states. First and final elements are always `0, T`. All elements lie in range `[0,T]`, see notes for behavior when field is set to an array with points outside this range.
- `precompute::TT`: Times for which the micromotion propagator and Floquet modes are precomputed. When initializing this struct, this field may be left blank, set to a Bool, or set to a list of timepoints. If set to `false`, no propagators will be stored except for the final period Hamiltonian. If left blank or set to `true`, this field will be set to the same value as `tlist`.  All elements  lie in range `[0,T]`. See notes for behavior when field is set to an array with points outside this range.
- `U_T::Qobj`: System propagator at time ``T``
- `Ulist::AbstractVector{Qobj}`: List of system propagators at times `tlist`
- `equasi::TE`: Time-independent quasienergies
"""
struct FloquetBasis{
    TT<:AbstractVector,
    TE<:AbstractVector,
    TQ<:AbstractVector{<:AbstractQuantumObject},
    TolT<:Real
}
    H::AbstractQuantumObject
    T::Real
    tlist::TT
    precompute::TT
    U_T::Qobj
    Ulist::TQ
    equasi::TE
    alg::AbstractODEAlgorithm
    abstol::TolT
    reltol::TolT

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
        alg::AbstractODEAlgorithm = Vern7(lazy = false),
        kwargs::Dict=Dict()
        ) where {TT<:AbstractVector}
        if T<=0
            throw(
                ArgumentError("`T` must be a nonzero positive real number")
            )
        end
        # enforce `tlist` and `precompute` rules
        tlist, precompute = _to_period_interval.([tlist, precompute], T) # enforce that all timepoints lie in interval [0,T)
        tlist = union(tlist, precompute) # ensure all times in precompute are in tlist
        tlist, precompute = [unique([0, tlist..., T]), unique([precompute..., T])] # ensure that period-propagator is calculated
        # solve for propagators
        kwargs[:saveat] = precompute
        sol = sesolve(H, qeye_like(H), tlist, alg=alg, kwargs=kwargs)
        sol_kwargs = sol.prob.kwargs
        Ulist = sol.states
        U_T = pop!(Ulist)
        # solve for quasienergies
        period_phases = eigenenergies(U_T)
        equasi = angle.(period_phases) ./ T

        new{typeof(tlist), typeof(equasi), typeof(Ulist)}(
            H,
            T,
            tlist,
            precompute,
            U_T,
            Ulist,
            equasi,
            sol.alg,
            sol_kwargs.abstol,
            sol_kwargs.reltol,
        )
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
function FloquetBasis(
    H::AbstractQuantumObject,
    T::Real,
    precompute::Bool=true;
    alg::AbstractODEAlgorithm = Vern7(lazy = false),
    kwargs::Dict=Dict()
    )
    tlist = range(0, T, 101)
    return FloquetBasis(H, T, tlist, precompute; alg=alg, kwargs=kwargs)
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
    tlist::AbstractVector,
    precompute::Bool=true;
    alg::AbstractODEAlgorithm = Vern7(lazy = false),
    kwargs::Dict=Dict())
    if precompute
        return FloquetBasis(H, T, tlist, tlist; alg=alg, kwargs=kwargs)
    else
        return FloquetBasis(H, T, tlist, Float64[]; alg=alg, kwargs=kwargs)
    end
end

function memoize_micromotion!(fb::FloquetBasis, tlist::AbstractVector, Ulist::AbstractVector{AbstractQuantumObject})
    length(tlist) == length(Ulist) ? nothing : throw(
        ArgumentError("tlist must be of same length as Ulist")
    )
    fb.precompute = [fb.precompute..., tlist...] |> unique |> sort
    insert_idxs = findall(x -> x∈tlist, fb.precompute)
    for (i_Ulist, i_fb) in enumerate(insert_idxs)
        insert!(fb.Ulist, i_fb, Ulist[i_Ulist])
    end
end

function propagator(fb::FloquetBasis, tlist::AbstractVector)
    # if tlist does not start at 0, propagate the system to tlist[1]
    #  at end of computation, the inverse of U(0, tlist[1]) will be right-multiplied to U(0, tlist[end])
    #  to obtain U(tlist[1], tlist[end])
    if tlist[1] != 0
        U_0 = propagator(fb, 0, tlist[1])
    else
        U_0 = qeye_like(fb.H)
    end
    kwargs = Dict([fb.kwargs... , :saveat=>:t_f])
    U_nT, Ulist_intra, _ = _compute_prop_list(fb, tlist, kwargs)
    return Ulist_intra[end] * U_nT * U_0^(-1)
end

function propagator!(fb::FloquetBasis, tlist::AbstractVector, kwargs::Dict=Dict())
    # tell _compute_prop_list to return micromotion propagators for every unmemoized timestep
    if tlist[1] !=0
        U_0 = propagator!(fb, 0, tlist[1], kwargs=kwargs)
    else
        U_0 = qeye_like(fb.H)
    end
    kwargs = Dict([kwargs..., :saveat=>:all_t])
    U_nT, Ulist_intra, tlist_intra = _compute_prop_list(fb, tlist, kwargs)
    memoize_micromotion!(fb, tlist_intra, Ulist_intra)
    return Ulist_intra[end] * U_nT * U_0^(-1)
end

function propagator(fb::FloquetBasis, t_0::Real, t_f::Real, n_timestep::Int64=100, kwargs::Dict=Dict())
    tlist = range(t_0, t_f, n_timestep)
    return propagator(fb, tlist, kwargs)
end

function propagator!(fb::FloquetBasis, t_0::Real, t_f::Real, n_timestep::Int64=100, kwargs::Dict=Dict())
    tlist = range(t_0, t_f, n_timestep)
    return propagator!(fb, tlist, kwargs)
end

function propagator(fb::FloquetBasis, t_f::Real, n_timestep::Int64=100, kwargs::Dict=Dict())
    return propagator(fb, 0, t_f, n_timestep, kwargs)
end

function propagator!(fb::FloquetBasis, t_f::Real, n_timestep::Int64=100, kwargs::Dict=Dict())
    return propagator!(fb, 0, t_f, n_timestep, kwargs)
end

function _compute_prop_list(fb::FloquetBasis, tlist::AbstractVector, kwargs::Dict)
    # find number of periods inside tlist
    nT, t_rem = fldmod(tlist[end], fb.T)
    # propagate by nT periods
    U_nT = fb.U_T^nT
    # propagate intra-period dynamics
    if t_rem==0 # if t_final is a multiple of T
        Ulist_intra = [qeye_like(fb.H)]
    elseif t_rem∈fb.precompute # if the correct micromotion propagator is memoized
        Ulist_intra = [fb.Ulist[findfirst(x->x==t_rem, fb.precompute)]]
    else # calculate micromotion operator from sesolve
        # calculate intra-period interval
        # TODO: Overide default progress bar to show number of micromotion propagators to memoize
        i_0 = findfirst(x -> x >= nT*fb.T, tlist)
        tlist_intra = _to_period_interval(tlist[i_0:end;], fb.T)
        tlist_intra =  tlist_intra[1] == 0 ? tlist_intra : [0, tlist_intra...]
        # Determine if all intra-period unitaries should be returned, or only final unitary
        kwargs[:saveat] = kwargs[:saveat] == :t_f ? (t_rem,) : tlist_intra
        Ulist_intra = sesolve(fb.H,
                              qeye_like(fb.H),
                              tlist,
                              alg=fb.alg,
                              kwargs=kwargs).states
    end
    return U_nT, Ulist_intra, tlist_intra
end


"""
    fsesolve(
    fb::FloquetBasis,
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector,
    alg::AbstractODEAlgorithm = Vern7(lazy=false),
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    params=NullParameters(),
    progress_bar::Union{Val, Bool} = Val(true),
    inplace::Union{Val,Bool}=Val(true),
    kwargs...;
    exact_t::Bool = false,
)

TBW
"""
function fsesolve(
    fb::FloquetBasis,
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector{TS},
    alg::AbstractODEAlgorithm = Vern7(lazy=false),
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    params=NullParameters(),
    progress_bar::Union{Val, Bool} = Val(true),
    inplace::Union{Val,Bool}=Val(true),
    kwargs...;
    efficient::Bool = false,
) where {TS<:Real}
    nsteps = length(tlist)
    if :saveat ∈ kwargs
        nstates = kwargs[:saveat]
    elseif isnothing(e_ops) || isempty(e_ops)
        kwargs[:saveat] = tlist
        nstates = length(tlist)
    else
        kwargs[:saveat] = tlist[end]
        nstates = 1
    end
    # pre-allocate solution memory
    sol = FloquetEvolutionSol(
        Vector{Float64}(undef, nsteps),
        Vector{Float64}(undef, nstates),
        Vector{QuantumObject{Ket}}(undef, nstates),
        isnothing(e_ops) || isempty(e_ops) ? nothing : Vector{Float64}(undef, nsteps),
        fb.alg,
        fb.abstol,
        fb.reltol
    )
    nothing
end
