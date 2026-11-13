export correlation_3op_2t, correlation_3op_1t, correlation_2op_2t, correlation_2op_1t

function _check_correlation_time_list(tlist::AbstractVector)
    any(t -> t == 0, tlist) ||
        throw(ArgumentError("The time list for calculating correlation function must contain the element `0`"))
    all(>=(0), tlist) ||
        throw(ArgumentError("All the elements in the time list for calculating correlation function must be positive."))
    issorted(tlist) || throw(ArgumentError("The time list for calculating correlation function must be sorted."))
    return nothing
end

@doc raw"""
    correlation_3op_2t(
        H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject,
        C::QuantumObject;
        alg::AbstractODEAlgorithm = DP5(),
        max_t_plus_tau::Union{Nothing,Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        matrix_form::Union{Val,Bool} = Val(false),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...
    )

Returns the two-time correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

# Notes

This function is implemented based on [`mesolve_map`](@ref) by the following adjustments:

- If the initial state `ψ0` is specified as `nothing`, then the [`steadystate`](@ref) will be used as the initial state.
- If `max_t_plus_tau` is specified, skip computation where ``t + \tau >`` `max_t_plus_tau`. Skipped entries are filled with ``0``. Default `nothing` means compute all entries (equivalent to `Inf`).
- The ``\tau``-evolution starting from every `t` in `tlist` is solved together as a single [`mesolve_map`](@ref), and `ensemblealg` (default to `EnsembleThreads()`) controls how that ensemble is parallelized.
- If `matrix_form = Val(true)`, the solving process is done with matrix form instead of vectorizing the density matrices, which is most useful for large Hilbert spaces, sparse operators, GPU calculations, and poorly sparse Liouvillians. See the `matrix_form` keyword argument of [`mesolve`](@ref) for more details.
"""
function correlation_3op_2t(
        H::AbstractQuantumObject{HOpType},
        ψ0::Union{Nothing, QuantumObject{StateOpType}},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        max_t_plus_tau::Union{Nothing, Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        matrix_form::Union{Val, Bool} = Val(false),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    # check tlist, τlist, and max_t_plus_tau
    _check_correlation_time_list(tlist)
    _check_correlation_time_list(τlist)
    max_time = something(max_t_plus_tau, Inf)
    max_time >= 0 || throw(ArgumentError("`max_t_plus_tau` must be `nothing` or a non-negative number."))

    matrix_form_val = makeVal(matrix_form)
    L = liouvillian(H, c_ops; matrix_form = matrix_form_val)
    if ψ0 isa Nothing
        ψ0 = steadystate(L)
    end

    # check the dimensions of Operator A, B, C, and ψ0
    # the dimension validity of L and ψ0 will be checked in mesolve
    check_mul_dimensions(A, ψ0)
    check_dimensions(A, B, C)

    # number of τ values with t + τ ≤ max_time for each t (same as QuTiP's searchsorted(..., side='right'))
    # tlist is sorted, so the t values without any τ left (t > max_time) are all at the end and are skipped
    nτ_list = [searchsortedlast(τlist, max_time - t) for t in tlist]
    tlist_kept = tlist[nτ_list .> 0]

    # generate OperatorKet for non-matrix form, so sol.states will store OperatorKet (skip redundant reshaping)
    ρ0 = getVal(matrix_form_val) ? ket2dm(ψ0) : mat2vec(ket2dm(ψ0))
    kwargs2 = merge((saveat = collect(tlist_kept),), (; kwargs...))

    # note that doesn't need to do the first mesolve if tlist_kept == [0] (for 2op_1t and 3op_1t cases)
    ρt_list =
        (tlist_kept == [0]) ? [ρ0] : mesolve(L, ρ0, tlist_kept; alg, progress_bar, matrix_form = matrix_form_val, kwargs2...).states

    # every t is one trajectory of a single mesolve_map ensemble, evolving C ρ(t) A over t + τ and recording ⟨B⟩
    # cache sprepost(C, A) to avoid repeated calls in the inner loop of mesolve_map
    C_A = cache_operator(sprepost(C, A; matrix_form = matrix_form_val), first(ρt_list))
    probs = [
        mesolveProblem(L, C_A * ρt_list[i], t .+ τlist[1:nτ_list[i]]; e_ops = [B], progress_bar = Val(false), matrix_form = matrix_form_val, kwargs...)
            for (i, t) in enumerate(tlist_kept)
    ]

    # each problem already has its own e_ops callback, so only the Liouvillian (its coefficients are updated in place) is copied
    _prob_func(_, ctx) = (prob = probs[ctx.sim_id].prob; remake(prob, f = _copy_for_trajectory(prob.f.f)))
    sols = mesolve_map(first(probs), probs, alg, ensemblealg; prob_func = _prob_func, safetycopy = false, progress_bar)

    # make the output correlation Matrix align with QuTiP
    # 1st dimension corresponds to tlist
    # 2nd dimension corresponds to τlist (the entries with t + τ > max_t_plus_tau are left as 0)
    corr = zeros(eltype(first(sols).expect), length(tlist), length(τlist))
    for (i, sol) in enumerate(sols)
        corr[i, 1:nτ_list[i]] .= @view sol.expect[1, :]
    end
    return corr
end

@doc raw"""
    correlation_3op_1t(
        H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject,
        C::QuantumObject;
        alg::AbstractODEAlgorithm = DP5(),
        matrix_form::Union{Val,Bool} = Val(false),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...
    )

Returns the two-time correlation function (with only one time coordinate ``\tau``) of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

# Notes

- If the initial state `ψ0` is specified as `nothing`, then the [`steadystate`](@ref) will be used as the initial state.
- If `matrix_form = Val(true)`, the solving process is done with matrix form instead of vectorizing the density matrices, which is most useful for large Hilbert spaces, sparse operators, GPU calculations, and poorly sparse Liouvillians. See the `matrix_form` keyword argument of [`mesolve`](@ref) for more details.
"""
function correlation_3op_1t(
        H::AbstractQuantumObject{HOpType},
        ψ0::Union{Nothing, QuantumObject{StateOpType}},
        τlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        matrix_form::Union{Val, Bool} = Val(false),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    # there is only a single t = 0, so the underlying mesolve_map ensemble has just one trajectory
    corr = correlation_3op_2t(H, ψ0, [0], τlist, c_ops, A, B, C; alg, ensemblealg = EnsembleSerial(), matrix_form = makeVal(matrix_form), progress_bar, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end

@doc raw"""
    correlation_2op_2t(
        H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject;
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool=false,
        max_t_plus_tau::Union{Nothing,Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        matrix_form::Union{Val,Bool} = Val(false),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...
    )

Returns the two-time correlation function of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

# Notes

- If the initial state `ψ0` is specified as `nothing`, then the [`steadystate`](@ref) will be used as the initial state.
- When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
- If `max_t_plus_tau` is specified, skip computation where ``t + \tau >`` `max_t_plus_tau`. Skipped entries are filled with ``0``. Default `nothing` means compute all entries (equivalent to `Inf`).
- The ``\tau``-evolution starting from every `t` in `tlist` is solved together as a single [`mesolve_map`](@ref), and `ensemblealg` (default to `EnsembleThreads()`) controls how that ensemble is parallelized. `alg` and `progress_bar` are forwarded to both the [`mesolve`](@ref) and [`mesolve_map`](@ref) calls.
- If `matrix_form = Val(true)`, the solving process is done with matrix form instead of vectorizing the density matrices, which is most useful for large Hilbert spaces, sparse operators, GPU calculations, and poorly sparse Liouvillians. See the `matrix_form` keyword argument of [`mesolve`](@ref) for more details.
"""
function correlation_2op_2t(
        H::AbstractQuantumObject{HOpType},
        ψ0::Union{Nothing, QuantumObject{StateOpType}},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool = false,
        max_t_plus_tau::Union{Nothing, Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        matrix_form::Union{Val, Bool} = Val(false),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    matrix_form_val = makeVal(matrix_form)
    C = eye(get_size(A.dimensions)[1], type = Operator(), dims = A.dims) # same as qeye_like(A), use A instead of H (cause H might be SuperOperator)
    if reverse
        corr = correlation_3op_2t(H, ψ0, tlist, τlist, c_ops, A, B, C; alg, max_t_plus_tau, ensemblealg, matrix_form = matrix_form_val, progress_bar, kwargs...)
    else
        corr = correlation_3op_2t(H, ψ0, tlist, τlist, c_ops, C, A, B; alg, max_t_plus_tau, ensemblealg, matrix_form = matrix_form_val, progress_bar, kwargs...)
    end

    return corr
end

@doc raw"""
    correlation_2op_1t(
        H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject;
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool=false,
        matrix_form::Union{Val,Bool} = Val(false),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...
    )

Returns the two-time correlation function (with only one time coordinate ``\tau``) of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

# Notes

- If the initial state `ψ0` is specified as `nothing`, then the [`steadystate`](@ref) will be used as the initial state.
- When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
- If `matrix_form = Val(true)`, the solving process is done with matrix form instead of vectorizing the density matrices, which is most useful for large Hilbert spaces, sparse operators, GPU calculations, and poorly sparse Liouvillians. See the `matrix_form` keyword argument of [`mesolve`](@ref) for more details.
"""
function correlation_2op_1t(
        H::AbstractQuantumObject{HOpType},
        ψ0::Union{Nothing, QuantumObject{StateOpType}},
        τlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool = false,
        matrix_form::Union{Val, Bool} = Val(false),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    # there is only a single t = 0, so the underlying mesolve_map ensemble has just one trajectory
    corr = correlation_2op_2t(H, ψ0, [0], τlist, c_ops, A, B; alg, reverse, ensemblealg = EnsembleSerial(), matrix_form = makeVal(matrix_form), progress_bar, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end
