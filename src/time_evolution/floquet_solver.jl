# script helper functions
function _to_period_interval(tlist::AbstractVector, T::Real)
    # function maps all elements ``t`` in `tlist` outside the interval ``[0, T)`` to an equivalent
    # time ``\tau`` such that ``mod(t, T) = \tau``
    if !isempty(tlist)
        tlist = mod.(tlist, T) |> unique |> sort
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
}
    times::TT1
    times_states::TT2
    states::TS
    expect::TE
    alg::AlgT
    abstol::TolT
    reltol::TolT
end


function _init_FloquetEvolutionSol(
    fb::FloquetBasis,
    tlist::AbstractVector{<:Real},
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    kwargs...
)
    # determine if expectation values need to be calculated
    has_eops = !(isnothing(e_ops) || isempty(e_ops))
    # determine timesteps at which to store propagated state
    nsteps = length(tlist)
    if haskey(kwargs, :saveat)
        # entry not checked by multiple dispatch, check if of correct type
        if !(kwargs[:saveat] isa AbstractVector{<:Real})
            throw(
                TypeError(:fsesolve, "saveat", AbstractVector{<:Real}, typeof(kwargs[:saveat]))
            )
        end
    elseif !has_eops
        kwargs[:saveat] = tlist
    else
        kwargs[:saveat] = [tlist[end]]
    end
    nstates = length(kwargs[:saveat])

    # determine whether expectation values must be collected
    #  if so, determine size
    n_eops = has_eops ? length(e_ops) : 0
    # pre-allocate solution memory
    sol = FloquetEvolutionSol(
        tlist,
        kwargs[:saveat],
        Vector{QuantumObject{Ket}}(undef, nstates),
        has_eops ? Array{ComplexF64}(undef, nsteps, n_eops) : nothing,
        fb.alg,
        fb.abstol,
        fb.reltol
    )
    return sol
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
- `Ulist::AbstractVector{Qobj}`: List of micromotion propagators at times `tlist`
- `equasi::TE`: Time-independent quasienergies
"""
struct FloquetBasis{
    TT<:AbstractVector{<:Real},
    TE<:AbstractVector{<:Real},
    TQ<:AbstractVector{<:AbstractQuantumObject},
    TolT<:Real
}
    H::AbstractQuantumObject
    T::Real
    precompute::TT
    U_T::Qobj
    Ulist::TQ
    equasi::TE
    alg::AbstractODEAlgorithm
    abstol::TolT
    reltol::TolT
    kwargs::Dict

    @doc raw"""
        FloquetBasis(H::AbstractQuantumObject, T::Real, tlist::AbstractVector{Real}, precompute::Bool = true; kwargs::Dict = Dict())

        DOCSTRING

        # Arguments:
        - `H`: Time-dependent system Hamiltonian.
        - `T`: Hamiltonian period such that ``\hat{H}(T+\tau_0) = \hat{H}(\tau_0)``.
        - `precompute`: Time vector containing points ``t`` at which to store the system propagator ``U(t)``.
        - `kwargs`: Additional keyword arguments to pass to ssesolve.

        # Notes:
        - If `precompute` contains elements outside the interval ``[0,T]``, a new time vector will be produced with all times ``t_k`` not in the interval mapped to an equivalent time ``\tau_k`` in the interval such that ``\hat{H}(t_k) = \hat{H}(\tau_k)``.

        # Returns:
        - `fbasis::FloquetBasis`: FloquetBasis object for the system evolving under the ``T``-periodic Hamiltonian `H`.
        """
    function FloquetBasis(
        H::AbstractQuantumObject,
        T::Real,
        precompute::TL;
        alg::AbstractODEAlgorithm = Vern7(lazy = false),
        kwargs::Dict=Dict()
        ) where {TL<:AbstractVector{<:Real}}
        if T<=0
            throw(
                ArgumentError("`T` must be a nonzero positive real number")
            )
        end
        # enforce `precompute` interval rule
        precompute = _to_period_interval(precompute, T)
        tlist = [0, precompute..., T] |> unique
        tlist = union(tlist, precompute) # ensure all times in precompute are in tlist
        # solve for propagators
        kwargs[:saveat] = precompute
        sol = sesolve(H, qeye_like(H), tlist, alg=alg, kwargs=kwargs)
        Ulist = sol.states
        U_T = pop!(Ulist)
        # solve for quasienergies
        period_phases = eigenenergies(U_T)
        equasi = angle.(period_phases) ./ T

        new{typeof(tlist), typeof(equasi), typeof(Ulist), typeof(sol.abstol)}(
            H,
            T,
            precompute,
            U_T,
            Ulist,
            equasi,
            sol.alg,
            sol.abstol,
            sol.reltol,
            kwargs,
        )
    end
end

@doc raw"""
    FloquetBasis(H::AbstractQuantumObject, T::Real; kwargs::Dict = Dict())

DOCSTRING

# Arguments:
- `H`: Time-dependent system Hamiltonian.
- `T`: Hamiltonian period such that ``\hat{H}(T+\tau_0) = \hat{H}(\tau_0)``.
- `kwargs`: Additional keyword arguments to pass to ssesolve.

# Notes
- Calling `FloquetBasis` without providing a time vector will create a FloquetBasis object with default `tlist=range(start:0, stop:T, length:101)`.

# Returns:
- `fbasis::FloquetBasis`: Floquet basis object for the system evolving under the time-dependent Hamiltonian `H`.
"""
function FloquetBasis(
    H::AbstractQuantumObject,
    T::TP,
    alg::AbstractODEAlgorithm = Vern7(lazy = false),
    kwargs::Dict=Dict()
    ) where {TP<:Real}
    precompute = Float64[]
    return FloquetBasis(H, T, precompute; alg=alg, kwargs=kwargs)
end


function memoize_micromotion!(
    fb::FloquetBasis,
    tlist::AbstractVector{<:Real};
    kwargs...
    )
    tlist = _to_period_interval(tlist, fb.T)
    propagator!(fb, tlist, kwargs...)
end

function memoize_micromotion!(fb::FloquetBasis, t::TP, U::QuantumObject{Operator}) where {TP<:Real}
    t = _to_period_interval(t, fb.T)
    t_idx = findfirst(x -> x>t, fb.precompute)
    t_idx = isnothing(t_idx) ? length(fb.precompute) + 1 : t_idx
    insert!(fb.precompute, t_idx, t)
    insert!(fb.Ulist, t_idx, U)
end



function memoize_micromotion!(
    fb::FloquetBasis,
    tlist::AbstractVector{<:Real},
    Ulist::AbstractVector{QuantumObject{Operator}}
    )
    for (t, U) in zip(tlist, Ulist)
        memoize_micromotion!(fb, t, U)
    end
end

function propagator(fb::FloquetBasis, t::TP; kwargs...) where {TP<:Real}
    return propagator(fb, 0, t, kwargs...)
end

function propagator!(fb::FloquetBasis, t::TP; kwargs...) where {TP<:Real}
    return propagator!(fb, 0, t; kwargs...)
end

function propagator(fb::FloquetBasis, t0::TP, tf::TP; kwargs...) where {TP<:Real}
    U0, U_nT, U_intra = _prop_list(fb, t0, tf; kwargs...)
    return U_intra * U_nT * U0'
end

function propagator!(fb::FloquetBasis, t0::TP, tf::TP; kwargs) where{TP<:Real}
    U0, U_nT, U_intra = _prop_list(fb, t0, tf; kwargs...)
    t0_T, tf_T = mod.([t0, tf], fb.T)
    for (tp, U) in [(t0_T, U0), (tf_T, U_intra)]
        if !(tp==0 || tp∈fb.precompute)
            memoize_micromotion!(fb, tp, U)
        end
    end
    return U_intra * U_nT * U0'
end

function _prop_list(fb::FloquetBasis, t0::TP, tf::TP; kwargs...) where {TP<:Real}
    U0 = (t0 == 0) ? qeye_like(fb.H) : propagator(fb, 0, t0, kwargs...)
    nT, t_rem = fldmod(tf, fb.T)
    U_nT = fb.U_T^nT
    if t_rem == 0
        U_intra = qeye_like(fb.H)
    elseif t_rem∈fb.precompute
        t_idx = findfirst(x->x==t_rem, fb.precompute)
        U_intra = fb.Ulist[t_idx]
    else
        if haskey(kwargs, :memoized_only) && kwargs[:memoized_only]
            throw(
                ErrorException(
                    "`memoized_only` keyword argument set to `true`, but "*
                        "the micromotion propagator corresponding to `tf=$tf` "*
                        "is not memoized in fb."
                )
            )
        end
        tlist = (0, t_rem)
        kwargs = Dict(fb.kwargs..., kwargs...)
        U_intra = sesolve(fb.H, qeye_like(fb.H), tlist; alg=fb.alg, kwargs...).states[end]
    end
    return U0, U_nT, U_intra
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
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    progress_bar::Union{Val, Bool} = Val(true),
    kwargs...,
) where {TS<:Real}
    return _fsesolve(fb, ψ0, tlist, propagator, e_ops, progress_bar; kwargs...)
end


function fsesolve!(
    fb::FloquetBasis,
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector{TS},
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    progress_bar::Union{Val, Bool} = Val(true),
    kwargs...,
) where {TS<:Real}
    return _fsesolve(fb, ψ0, tlist, propagator!, e_ops, progress_bar; kwargs...)
end


function _fsesolve(
    fb::FloquetBasis,
    ψ0::QuantumObject{Ket},
    tlist::AbstractVector{TS},
    pfunc::Function,
    e_ops::Union{Nothing, AbstractVector, Tuple} = nothing,
    progress_bar::Union{Val, Bool} = Val(true),
    kwargs...,
) where {TS<:Real}
    sol = _init_FloquetEvolutionSol(
        fb,
        tlist,
        e_ops,
        kwargs...
    )
    # convert progress_bar to boolean type
    if !(progress_bar isa Bool)
        progress_bar = typeof(progress_bar).parameters[1]
    end
    if progress_bar
        pbar = Progress(length(tlist), showspeed=true)
    end
    for (step, t) in enumerate(tlist)
        if t==0
            U = qeye_like(fb.H)
        else
            U = pfunc(fb, 0, t; kwargs...)
        end
        ψt = U * ψ0
        if haskey(kwargs, :saveat)
            state_idx = findfirst(x->x==t, kwargs[:saveat])
            if !isnothing(state_idx)
                sol.states[state_idx] = ψt
            end
        end
        if !isnothing(sol.expect)
            sol.expect[step,:] = expect.(e_ops, ψt)
        end
        progress_bar ? next!(pbar) : nothing
    end
    progress_bar ? finish!(pbar) : nothing
    return sol
end

function _state_list(fb::FloquetBasis)
    phases = exp.( fb.T * 1im .* fb.equasi) |>
        Diagonal |>
        x -> Qobj(x, dims=fb.U_T.dims)
    U = phases * fb.U_T
    _, ψ0, U0 = eigenstates(U)
    return ψ0, U0.data
end

function _state_list(fb::FloquetBasis, t::TP, pfunc::Function) where {TP<:Real}
    _, ψ0, U0_mat = _mode_list(fb)
    Ut_mat = (pfunc(fb, t).data * U0_mat)
    Ut_list = eachcol(Ut_mat) |> collect
    ψt = Qobj.(Ut_list, dims=ψ0[1].dims)
    return ψt, Ut_mat
end

"""
TODO: mode, state, from_floquet_basis, to_floquet_basis
"""
