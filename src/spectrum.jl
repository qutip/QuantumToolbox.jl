export spectrum, spectrum_correlation_fft
export SpectrumSolver, ExponentialSeries, PseudoInverse

abstract type SpectrumSolver end

@doc raw"""
    ExponentialSeries(; tol = 1e-14, calc_steadystate = false)

A solver which solves [`spectrum`](@ref) by finding the eigen decomposition of the Liouvillian [`SuperOperator`](@ref) and calculate the exponential series.
"""
struct ExponentialSeries{T<:Real,CALC_SS} <: SpectrumSolver
    tol::T
    ExponentialSeries(tol::T, calc_steadystate::Bool = false) where {T} = new{T,calc_steadystate}(tol)
end

ExponentialSeries(; tol = 1e-14, calc_steadystate = false) = ExponentialSeries(tol, calc_steadystate)

@doc raw"""
    PseudoInverse(; alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES())

A solver which solves [`spectrum`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the `alg`orithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/).
"""
struct PseudoInverse{MT<:SciMLLinearSolveAlgorithm} <: SpectrumSolver
    alg::MT
end

PseudoInverse(; alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES()) = PseudoInverse(alg)

@doc raw"""
    spectrum(H::QuantumObject,
        ωlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
        B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject};
        solver::SpectrumSolver=ExponentialSeries(),
        kwargs...)

Calculate the spectrum of the correlation function

```math
S(\omega) = \int_{-\infty}^\infty \lim_{t \rightarrow \infty} \left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle e^{-i \omega \tau} d \tau
```

See also the following list for `SpectrumSolver` docstrings:
- [`ExponentialSeries`](@ref)
- [`PseudoInverse`](@ref)
"""
function spectrum(
    H::QuantumObject{MT1,HOpType},
    ωlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject};
    solver::SpectrumSolver = ExponentialSeries(),
    kwargs...,
) where {MT1<:AbstractMatrix,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    return _spectrum(liouvillian(H, c_ops), ωlist, A, B, solver; kwargs...)
end

function _spectrum_get_rates_vecs_ss(L, solver::ExponentialSeries{T,true}) where {T}
    result = eigen(L)
    rates, vecs = result.values, result.vectors

    return rates, vecs, steadystate(L).data
end

function _spectrum_get_rates_vecs_ss(L, solver::ExponentialSeries{T,false}) where {T}
    result = eigen(L)
    rates, vecs = result.values, result.vectors

    ss_idx = findmin(abs2, rates)[2]
    ρss = vec2mat(@view(vecs[:, ss_idx]))
    ρss = (ρss + ρss') / 2
    ρss ./= tr(ρss)

    return rates, vecs, ρss
end

function _spectrum(
    L::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    ωlist::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    solver::ExponentialSeries;
    kwargs...,
) where {T1,T2,T3}
    allequal((L.dims, A.dims, B.dims)) ||
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))

    rates, vecs, ρss = _spectrum_get_rates_vecs_ss(L, solver)

    ρ0 = B.data * ρss
    v = vecs \ mat2vec(ρ0)

    amps = map(i -> v[i] * tr(A.data * vec2mat(@view(vecs[:, i]))), eachindex(rates))
    idxs = findall(x -> abs(x) > solver.tol, amps)
    amps, rates = amps[idxs], rates[idxs]

    # spec = map(ω -> 2 * real(sum(@. amps * (1 / (1im * ω - rates)))), ωlist)
    amps_rates = zip(amps, rates)
    spec = map(ω -> 2 * real(sum(x -> x[1] / (1im * ω - x[2]), amps_rates)), ωlist)

    return spec
end

function _spectrum(
    L::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    ωlist::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    solver::PseudoInverse;
    kwargs...,
) where {T1,T2,T3}
    allequal((L.dims, A.dims, B.dims)) ||
        throw(DimensionMismatch("The quantum objects are not of the same Hilbert dimension."))

    ωList = convert(Vector{_FType(L)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    spec = Vector{_FType(L)}(undef, Length)

    # calculate vectorized steadystate, multiply by operator B on the left (spre)
    ρss = mat2vec(steadystate(L))
    b = (spre(B) * ρss).data

    # multiply by operator A on the left (spre) and then perform trace operation
    D = prod(L.dims)
    _tr = SparseVector(D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(_CType(L), D)) # same as vec(system_identity_matrix)
    _tr_A = transpose(_tr) * spre(A).data

    I_cache = I(D^2)

    ω = popfirst!(ωList)
    cache = init(LinearProblem(L.data - 1im * ω * I_cache, b), solver.alg, kwargs...)
    sol = solve!(cache)
    spec[1] = -2 * real(dot(_tr_A, sol.u))

    for (idx, ω) in enumerate(ωList)
        cache.A = L.data - 1im * ω * I_cache
        sol = solve!(cache)

        # trace over the Hilbert space of system (expectation value)
        spec[idx+1] = -2 * real(dot(_tr_A, sol.u))
    end

    return spec
end

@doc raw"""
    spectrum_correlation_fft(tlist, corr; inverse=Val(false))

Calculate the power spectrum corresponding to a two-time correlation function using fast Fourier transform (FFT).

# Parameters
- `tlist::AbstractVector`: List of times at which the two-time correlation function is given.
- `corr::AbstractVector`: List of two-time correlations corresponding to the given time point in `tlist`.
- `inverse::Union{Val,Bool}`: Whether to use the inverse Fourier transform or not. Default to `Val(false)`.

# Returns
- `ωlist`: the list of angular frequencies ``\omega``.
- `Slist`: the list of the power spectrum corresponding to the angular frequencies in `ωlist`.
"""
function spectrum_correlation_fft(tlist::AbstractVector, corr::AbstractVector; inverse::Union{Val,Bool} = Val(false))
    N = length(tlist)
    dt_list = diff(tlist)
    dt = dt_list[1]

    all(≈(dt), dt_list) || throw(ArgumentError("tlist must be equally spaced for FFT."))

    # power spectrum list
    F = getVal(inverse) ? N * ifft(corr) : fft(corr)
    Slist = 2 * dt * real(fftshift(F))

    # angular frequency list
    ωlist = 2 * π * fftshift(fftfreq(N, 1 / dt))

    return ωlist, Slist
end
