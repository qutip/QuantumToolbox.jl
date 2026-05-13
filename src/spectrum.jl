export spectrum, spectrum_correlation_fft
export SpectrumSolver, ExponentialSeries, PseudoInverse, Lanczos

abstract type SpectrumSolver end

@doc raw"""
    ExponentialSeries(; tol = 1e-14, calc_steadystate = false)

A solver which solves [`spectrum`](@ref) by finding the eigen decomposition of the Liouvillian [`SuperOperator`](@ref) and calculate the exponential series.
"""
struct ExponentialSeries{T <: Real, CALC_SS} <: SpectrumSolver
    tol::T
    ExponentialSeries(tol::T, calc_steadystate::Bool = false) where {T} = new{T, calc_steadystate}(tol)
end

ExponentialSeries(; tol = 1.0e-14, calc_steadystate = false) = ExponentialSeries(tol, calc_steadystate)

@doc raw"""
    PseudoInverse(; alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES())

A solver which solves [`spectrum`](@ref) by finding the inverse of Liouvillian [`SuperOperator`](@ref) using the `alg`orithms given in [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/).
"""
struct PseudoInverse{MT <: SciMLLinearSolveAlgorithm} <: SpectrumSolver
    alg::MT
end

PseudoInverse(; alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES()) = PseudoInverse(alg)

@doc raw"""
    Lanczos(; tol = 1e-8, maxiter = 5000, verbose = 0)

A solver which solves [`spectrum`](@ref) by using a non-symmetric Lanczos variant of the algorithm in [Koch2011](https://www.cond-mat.de/events/correl11/manuscripts/koch.pdf).
The nonsymmetric Lanczos algorithm is adapted from Algorithm 6.6 in [Saad2011](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf).
The running estimate is updated via a [Wallis-Euler recursion](https://en.wikipedia.org/wiki/Continued_fraction).
"""
Base.@kwdef struct Lanczos{T <: Real, SS <: Union{Nothing, <:SteadyStateSolver}} <: SpectrumSolver
    tol::T = 1.0e-8
    maxiter::Int = 5000
    verbose::Int = 0
    steadystate_solver::SS = nothing
end

@doc raw"""
    spectrum(H::QuantumObject,
        ωlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        solver::SpectrumSolver=ExponentialSeries(),
        kwargs...)

Calculate the spectrum of the correlation function

```math
S(\omega) = \int_{-\infty}^\infty \lim_{t \rightarrow \infty} \left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle e^{-i \omega \tau} d \tau
```

See also the following list for `SpectrumSolver` docstrings:
- [`ExponentialSeries`](@ref)
- [`PseudoInverse`](@ref)
- [`Lanczos`](@ref)
"""
function spectrum(
        H::QuantumObject{HOpType},
        ωlist::AbstractVector,
        c_ops::Union{Nothing, AbstractVector, Tuple},
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        solver::SpectrumSolver = ExponentialSeries(),
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}}

    !isendomorphic(H.dimensions) && _non_endomorphic_dims_error("Hamiltonian or Liouvillian for spectrum", H.dimensions)

    L = liouvillian(H, c_ops)
    (L isa QuantumObject) || throw(ArgumentError("spectrum only supports (time-independent) QuantumObject in c_ops"))
    check_mul_dimensions(L, A)
    check_dimensions(A, B)

    return _spectrum(L, ωlist, A, B, solver; kwargs...)
end

function _spectrum_get_rates_vecs_ss(L, solver::ExponentialSeries{T, true}) where {T}
    result = eigen(L)
    rates, vecs = result.values, result.vectors

    return rates, vecs, steadystate(L).data
end

function _spectrum_get_rates_vecs_ss(L, solver::ExponentialSeries{T, false}) where {T}
    result = eigen(L)
    rates, vecs = result.values, result.vectors

    ss_idx = findmin(abs2, rates)[2]
    ρss = vec2mat(@view(vecs[:, ss_idx]))
    ρss = (ρss + ρss') / 2
    ρss ./= tr(ρss)

    return rates, vecs, ρss
end

function _spectrum(
        L::QuantumObject{SuperOperator},
        ωlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        solver::ExponentialSeries;
        kwargs...,
    )
    rates, vecs, ρss = _spectrum_get_rates_vecs_ss(L, solver)

    ρ0 = B.data * ρss
    v = vecs \ mat2vec(ρ0)

    amps = map(i -> v[i] * tr(A.data * vec2mat(@view(vecs[:, i]))), eachindex(rates))
    idxs = findall(x -> abs(x) > solver.tol, amps)
    amps, rates = amps[idxs], rates[idxs]

    amps_rates = zip(amps, complex.(rates))
    spec = map(ω -> 2 * real(sum(x -> x[1] / (1im * ω - x[2]), amps_rates)), ωlist)

    return spec
end

function _spectrum(
        L::QuantumObject{SuperOperator},
        ωlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        solver::PseudoInverse;
        kwargs...,
    )
    ωList = convert(Vector{_float_type(L)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    spec = Vector{_float_type(L)}(undef, Length)

    # calculate vectorized steadystate, multiply by operator B on the left (spre)
    ρss = mat2vec(steadystate(L))
    b = (spre(B) * ρss).data

    # multiply by operator A on the left (spre) and then perform trace operation
    D = size(A, 1)
    _tr = SparseVector(D^2, [1 + n * (D + 1) for n in 0:(D - 1)], ones(_complex_float_type(L), D)) # same as vec(system_identity_matrix)
    _tr_A = transpose(_tr) * spre(A).data

    Id = Eye(D^2)

    # DO the idx = 1 case
    ω = ωList[1]
    cache = init(LinearProblem(L.data - 1im * ω * Id, b), solver.alg, kwargs...)
    sol = solve!(cache)
    spec[1] = -2 * real(dot(_tr_A, sol.u))
    popfirst!(ωList)
    for (idx, ω) in enumerate(ωList)
        cache.A = L.data - 1im * ω * Id
        sol = solve!(cache)

        # trace over the Hilbert space of system (expectation value)
        spec[idx + 1] = -2 * real(dot(_tr_A, sol.u))
    end

    return spec
end

function _spectrum(
        L::QuantumObject{SuperOperator},
        ωlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        solver::Lanczos,
    )
    # Define type shortcuts
    fT = _float_type(L)
    cT = _complex_float_type(L)

    # Calculate |v₁> = B|ρss>
    ρss =
        isnothing(solver.steadystate_solver) ? mat2vec(steadystate(L)) :
        mat2vec(steadystate(L; solver = solver.steadystate_solver))
    vₖ = (spre(B) * ρss).data

    # Define (possibly GPU) vector type
    vT = typeof(vₖ)

    # Calculate <w₁| = <I|A
    D = size(A, 1)
    Ivec = SparseVector(D^2, [1 + n * (D + 1) for n in 0:(D - 1)], ones(cT, D)) # same as vec(system_identity_matrix)
    wₖ = spre(A).data' * vT(Ivec)

    # Store the normalization factor for the Green's function before renormalizing |v₁> and <w₁|
    gfNorm = dot(wₖ, vₖ)
    if gfNorm ≈ 0.0
        throw(AssertionError("⟨w₀|v₀⟩ = 0, please check your A and B operators."))
    end
    scalingF = sqrt(abs(gfNorm))
    vₖ ./= scalingF
    wₖ ./= conj(gfNorm / scalingF)

    # Handle input frequency range
    ωList = vT(convert(Vector{fT}, ωlist))  # Make sure they're real frequencies and potentially on GPU
    Length = length(ωList)

    # Current and previous estimates of the spectrum
    lanczosFactor = vT(zeros(cT, Length))
    lanczosFactor₋₁ = vT(zeros(cT, Length))

    # Tridiagonal matrix elements
    αₖ = cT(0)
    βₖ = cT(-1)
    δₖ = fT(+1)
    βₖδₖ = βₖ * δₖ

    # Current and up to second-to-last A and B Euler sequences
    A₋₂ = vT(ones(cT, Length))
    A₋₁ = vT(zeros(cT, Length))
    Aₖ = vT(zeros(cT, Length))
    B₋₂ = vT(zeros(cT, Length))
    B₋₁ = vT(ones(cT, Length))
    Bₖ = vT(zeros(cT, Length))

    # Maximum norm and residue
    maxNorm = vT(zeros(cT, length(ωList)))
    maxResidue = fT(0.0)

    # Previous and next left/right Krylov vectors
    v₋₁ = vT(zeros(cT, D^2))
    v₊₁ = vT(zeros(cT, D^2))
    w₋₁ = vT(zeros(cT, D^2))
    w₊₁ = vT(zeros(cT, D^2))

    # Frequency of renormalization
    renormFrequency = 1

    # Loop over the Krylov subspace(s)
    for k in 1:solver.maxiter
        # k-th diagonal element
        mul!(w₊₁, L.data', wₖ)
        αₖ = dot(w₊₁, vₖ)

        # Update A(k), B(k) and continuous fraction; normalization avoids overflow
        Aₖ .= (-1im .* ωList .+ αₖ) .* A₋₁ .- βₖδₖ .* A₋₂
        Bₖ .= (-1im .* ωList .+ αₖ) .* B₋₁ .- βₖδₖ .* B₋₂
        lanczosFactor₋₁ .= lanczosFactor
        lanczosFactor .= Aₖ ./ Bₖ

        # Renormalize Euler sequences to avoid overflow
        if k % renormFrequency == 0
            map!((x, y) -> max(abs(x), abs(y)), maxNorm, Aₖ, Bₖ)
            Aₖ ./= maxNorm
            Bₖ ./= maxNorm
            A₋₁ ./= maxNorm
            B₋₁ ./= maxNorm
        end

        # Check for convergence

        residueNorm = max(maximum(abs, lanczosFactor), maximum(abs, lanczosFactor₋₁))
        lanczosFactor₋₁ .-= lanczosFactor
        maxResidue = maximum(abs, lanczosFactor₋₁) / residueNorm
        if maxResidue <= solver.tol
            if solver.verbose > 1
                println("spectrum(): solver::Lanczos converged after $(k) iterations")
            end
            break
        end

        # (k+1)-th left/right vectors, orthogonal to previous ones
        mul!(v₊₁, L.data, vₖ)
        v₊₁ .= v₊₁ .- αₖ .* vₖ .- βₖ .* v₋₁
        w₊₁ .= w₊₁ .- conj(αₖ) .* wₖ .- conj(δₖ) .* w₋₁
        v₋₁, vₖ = vₖ, v₋₁
        vₖ, v₊₁ = v₊₁, vₖ
        w₋₁, wₖ = wₖ, w₋₁
        wₖ, w₊₁ = w₊₁, wₖ

        # k-th off-diagonal elements
        βₖδₖ = dot(wₖ, vₖ)
        if βₖδₖ ≈ 0.0
            if solver.verbose > 0
                @warn "spectrum(): solver::Lanczos experienced orthogonality breakdown after $(k) iterations"
                @warn "spectrum(): βₖδₖ = $(βₖδₖ)"
            end
            break
        end
        δₖ = sqrt(abs(βₖδₖ))
        βₖ = βₖδₖ / δₖ

        # Normalize (k+1)-th left/right vectors
        vₖ ./= δₖ
        wₖ ./= conj(βₖ)

        # Update everything for the next cycle
        A₋₂, A₋₁ = A₋₁, A₋₂
        A₋₁, Aₖ = Aₖ, A₋₁
        B₋₂, B₋₁ = B₋₁, B₋₂
        B₋₁, Bₖ = Bₖ, B₋₁
    end

    if solver.verbose > 0 && maxResidue > solver.tol
        println("spectrum(): maxiter = $(solver.maxiter) reached before convergence!")
        println("spectrum(): Max residue = $maxResidue")
        println("spectrum(): Consider increasing maxiter and/or tol")
    end

    # Restore the norm
    lanczosFactor .= gfNorm .* lanczosFactor

    return -2 .* real(lanczosFactor)
end

@doc raw"""
    spectrum_correlation_fft(tlist, corr; inverse=false)

Calculate the power spectrum corresponding to a two-time correlation function using fast Fourier transform (FFT).

# Parameters
- `tlist::AbstractVector`: List of time points at which the two-time correlation function is given.
- `corr::AbstractVector`: List of two-time correlations corresponding to the given time point in `tlist`.
- `inverse::Bool`: Whether to use the inverse Fourier transform or not. Default to `false`.

# Returns
- `ωlist`: the list of angular frequencies ``\omega``.
- `Slist`: the list of the power spectrum corresponding to the angular frequencies in `ωlist`.
"""
function spectrum_correlation_fft(tlist::AbstractVector, corr::AbstractVector; inverse::Bool = false)
    N = length(tlist)
    dt_list = diff(tlist)
    dt = dt_list[1]

    all(≈(dt), dt_list) || throw(ArgumentError("tlist must be equally spaced for FFT."))

    # power spectrum list
    F = inverse ? N * ifft(corr) : fft(corr)
    Slist = 2 * dt * real(fftshift(F))

    # angular frequency list
    ωlist = 2 * π * fftshift(fftfreq(N, 1 / dt))

    return ωlist, Slist
end
