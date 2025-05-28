export spectrum, spectrum_correlation_fft
export SpectrumSolver, ExponentialSeries, PseudoInverse, Lanczos

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
    Lanczos(; tol = 1e-8, maxiter = 5000, verbose = 0)

A solver which solves [`spectrum`](@ref) by using a non-symmetric Lanczos variant of the algorithm in [Koch2011](https://www.cond-mat.de/events/correl11/manuscripts/koch.pdf).
The nonsymmetric Lanczos algorithm is adapted from Algorithm 6.6 in [Saad2011](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf).
The running estimate is updated via a [Wallis-Euler recursion](https://en.wikipedia.org/wiki/Continued_fraction).
"""
struct Lanczos{T<:Real,IT<:Int} <: SpectrumSolver
    tol::T
    maxiter::IT
    verbose::IT
end

Lanczos(; tol = 1e-8, maxiter = 5000, verbose = 0) = Lanczos(tol, maxiter, verbose)

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
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{Operator},
    B::QuantumObject{Operator};
    solver::SpectrumSolver = ExponentialSeries(),
    kwargs...,
) where {HOpType<:Union{Operator,SuperOperator}}
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
    L::QuantumObject{SuperOperator},
    ωlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    solver::ExponentialSeries;
    kwargs...,
)
    check_dimensions(L, A, B)

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
    L::QuantumObject{SuperOperator},
    ωlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    solver::PseudoInverse;
    kwargs...,
)
    check_dimensions(L, A, B)

    ωList = convert(Vector{_FType(L)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    spec = Vector{_FType(L)}(undef, Length)

    # calculate vectorized steadystate, multiply by operator B on the left (spre)
    ρss = mat2vec(steadystate(L))
    b = (spre(B) * ρss).data

    # multiply by operator A on the left (spre) and then perform trace operation
    D = prod(L.dimensions)
    _tr = SparseVector(D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(_CType(L), D)) # same as vec(system_identity_matrix)
    _tr_A = transpose(_tr) * spre(A).data

    Id = I(D^2)

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
        spec[idx+1] = -2 * real(dot(_tr_A, sol.u))
    end

    return spec
end

function _spectrum(
    L::QuantumObject{SuperOperator},
    ωlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    solver::Lanczos;
    kwargs...,
)
    check_dimensions(L, A, B)

    # Define type shortcuts
    fT = _FType(L)
    cT = _CType(L)

    # Handle input frequency range
    ωList = convert(Vector{fT}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    #spec = Vector{fT}(undef, Length)

    # Calculate |v₁> = B|ρss>
    ρss = mat2vec(steadystate(L))
    vₖ = Array{cT}((spre(B) * ρss).data)

    # Calculate <w₁| = <I|A
    D = prod(L.dimensions)
    Ivec = SparseVector(D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(cT, D)) # same as vec(system_identity_matrix)
    wₖ = transpose(typeof(vₖ)(Ivec)) * spre(A).data

    # Store the norm of the Green's function before renormalizing |v₁> and <w₁|
    gfNorm = abs(wₖ * vₖ)
    vₖ ./= sqrt(gfNorm)
    wₖ ./= sqrt(gfNorm)

    # println("  type: $(typeof(vₖ))")
    # println("  type: $(typeof(wₖ))")

    # Current and previous estimates of the spectrum
    lanczosFactor   = zeros(cT, Length)
    lanczosFactor₋₁ = zeros(cT, Length)

    # Tridiagonal matrix elements
    αₖ = cT( 0)
    βₖ = cT(-1)
    δₖ = cT(+1)

    # Current and up to second-to-last A and B Euler sequences
    A₋₂ =  ones(cT, Length)
    A₋₁ = zeros(cT, Length)
    Aₖ  = zeros(cT, Length)
    B₋₂ = zeros(cT, Length)
    B₋₁ =  ones(cT, Length)
    Bₖ  = zeros(cT, Length)

    # Maximum norm and residue
    maxNorm    = zeros(cT, length(ωList))
    maxResidue = fT(0.0)

    # Previous and next left/right Krylov vectors
    v₋₁ = zeros(cT, (D^2, 1))
    v₊₁ = zeros(cT, (D^2, 1))
    w₋₁ = zeros(cT, (1, D^2))
    w₊₁ = zeros(cT, (1, D^2))

    # Frequency of renormalization
    renormFrequency::typeof(solver.maxiter) = 1

    # Loop over the Krylov subspace(s)
    for k in 1:solver.maxiter
        # k-th diagonal element
        w₊₁ = wₖ * L.data
        αₖ = w₊₁ * vₖ
        
        # Update A(k), B(k) and continuous fraction; normalization avoids overflow
        Aₖ .= (-1im .* ωList .+ αₖ) .* A₋₁ .- (βₖ * δₖ) .* A₋₂
        Bₖ .= (-1im .* ωList .+ αₖ) .* B₋₁ .- (βₖ * δₖ) .* B₋₂
        lanczosFactor₋₁ .= lanczosFactor
        lanczosFactor .= Aₖ ./ Bₖ

        # Renormalize Euler sequences to avoid overflow
        if k % renormFrequency == 0
            maxNorm .= max.(abs.(Aₖ), abs.(Bₖ))  # Note: the MATLAB and C++ codes return the actual complex number
            Aₖ ./= maxNorm
            Bₖ ./= maxNorm
            A₋₁ ./= maxNorm
            B₋₁ ./= maxNorm
        end

        # Check for convergence
        maxResidue = maximum(abs.(lanczosFactor .- lanczosFactor₋₁)) / 
                     max(maximum(abs.(lanczosFactor)), maximum(abs.(lanczosFactor₋₁)))
        if maxResidue <= solver.tol
            if solver.verbose > 1
                println("spectrum(): solver::Lanczos converged after $(k) iterations")
            end
            break
        end

        # (k+1)-th left/right vectors, orthogonal to previous ones
        # Consider using explicit BLAS calls
        v₊₁ = L.data * vₖ
        v₊₁ .= v₊₁ .- αₖ .* vₖ .- βₖ .* v₋₁
        w₊₁ .= w₊₁ .- αₖ .* wₖ .- δₖ .* w₋₁
        v₋₁ .= vₖ
        w₋₁ .= wₖ
        vₖ  .= v₊₁
        wₖ  .= w₊₁

        # k-th off-diagonal elements
        buf = wₖ * vₖ
        δₖ  = sqrt(abs(buf))
        βₖ  = buf / δₖ

        # Normalize (k+1)-th left/right vectors
        vₖ ./= δₖ
        wₖ ./= βₖ

        # Update everything for the next cycle
        A₋₂ .= A₋₁
        A₋₁ .= Aₖ
        B₋₂ .= B₋₁
        B₋₁ .= Bₖ
    end

    if solver.verbose > 0
        if maxResidue > solver.tol
            println("spectrum(): maxiter = $(solver.maxiter) reached before convergence!")
            println("spectrum(): Max residue = $maxResidue")
            println("spectrum(): Consider increasing maxiter and/or tol")
        end
    end

    # Restore the norm
    lanczosFactor .= gfNorm .* lanczosFactor

    return -2 .* real( lanczosFactor )
end

@doc raw"""
    spectrum_correlation_fft(tlist, corr; inverse=false)

Calculate the power spectrum corresponding to a two-time correlation function using fast Fourier transform (FFT).

# Parameters
- `tlist::AbstractVector`: List of times at which the two-time correlation function is given.
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
