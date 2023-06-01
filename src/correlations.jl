abstract type SpectrumSolver end
struct FFTCorrelation <: SpectrumSolver end
struct ExponentialSeries <: SpectrumSolver
    tol::Real
end

ExponentialSeries(;tol=1e-14) = ExponentialSeries(tol)

@doc raw"""
    correlation_3op_2t(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector,
        τ_l::AbstractVector,
        A::QuantumObject,
        B::QuantumObject,
        C::QuantumObject,
        c_ops::AbstractVector=[];
        kwargs...)

Returns the two-times correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``:
``\expval{\hat{A}(t) \hat{B}(t + \tau) \hat{C}(t)}`` for a given initial state ``\ket{\psi_0}``.
"""
function correlation_3op_2t(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    C::QuantumObject{<:AbstractArray{T5},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    kwargs...) where {T1,T2,T3,T4,T5,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    (H.dims == ψ0.dims && H.dims == A.dims &&
    H.dims == B.dims && H.dims == C.dims) || throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    kwargs2 = merge(kwargs, Dict(:saveat => collect(t_l)))
    ρt = mesolve(H, ψ0, t_l, c_ops; kwargs2...).states

    corr = map((t,ρ)->mesolve(H, C*ρ*A, τ_l .+ t, c_ops, e_ops=[B]; kwargs...).expect[1,:], t_l, ρt)

    corr
end

@doc raw"""
    correlation_2op_2t(H::QuantumObject,
        ψ0::QuantumObject,
        t_l::AbstractVector,
        τ_l::AbstractVector,
        A::QuantumObject,
        B::QuantumObject,
        c_ops::AbstractVector=[];
        reverse::Bool=false,
        kwargs...)

Returns the two-times correlation function of two operators ``\hat{A}`` and ``\hat{B}``` 
at different times ``\expval{\hat{A}(t + \tau) \hat{B}(t)}``.
When ``reverse=true``, the correlation function is calculated as ``\expval{\hat{A}(t) \hat{B}(t + \tau)}``.
"""
function correlation_2op_2t(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector,
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    reverse::Bool=false,
    kwargs...) where {T1,T2,T3,T4,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    C = eye(prod(H.dims), dims=H.dims)
    if reverse
        corr = correlation_3op_2t(H, ψ0, t_l, τ_l, A, B, C, c_ops; kwargs...)
    else
        corr = correlation_3op_2t(H, ψ0, t_l, τ_l, C, A, B, c_ops; kwargs...)
    end
    
    reduce(hcat, corr)
end

@doc raw"""
    correlation_2op_1t(H::QuantumObject,
        ψ0::QuantumObject,
        τ_l::AbstractVector,
        A::QuantumObject,
        B::QuantumObject,
        c_ops::AbstractVector=[];
        reverse::Bool=false,
        kwargs...)

Returns the one-time correlation function of two operators ``\hat{A}`` and ``\hat{B}`` 
at different times ``\expval{\hat{A}(\tau) \hat{B}(0)}``.
When ``reverse=true``, the correlation function is calculated as ``\expval{\hat{A}(0) \hat{B}(\tau)}``.
"""
function correlation_2op_1t(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    reverse::Bool=false,
    kwargs...) where {T1,T2,T3,T4,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    corr = correlation_2op_2t(H, ψ0, [0], τ_l, A, B, c_ops; reverse=reverse, kwargs...)
    
    corr[:,1]
end

@doc raw"""
    spectrum(H::QuantumObject,
        ω_max::Real, Nsamples::Integer,
        A::QuantumObject,
        B::QuantumObject,
        c_ops::AbstractVector=[];
        alg=Vern7(),
        H_t=nothing,
        params::AbstractVector=[],
        progress::Bool=true,
        callbacks=[],
        kwargs...)

Returns the emission spectrum ``S(\omega) = \int_{-\infty}^\infty \expval{\hat{A}(\tau) \hat{B}(0)} e^{-i \omega \tau} d \tau``.
"""
function spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_max::Real, Nsamples::Integer,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    solver::MySolver=ExponentialSeries(),
    kwargs...) where {T1,T2,T3,
            HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}, MySolver<:SpectrumSolver}
    
    return _spectrum(H, ω_max, Nsamples, A, B, c_ops, solver; kwargs...)
end

function _spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_max::Real, Nsamples::Integer,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector,
    solver::FFTCorrelation;
    kwargs...) where {T1,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    
    ω_max = abs(ω_max)
    dω = 2*ω_max/(Nsamples-1)
    ω_l = -ω_max:dω:ω_max

    T = 2π/(ω_l[2]-ω_l[1])
    τ_l = range(0,T, length=length(ω_l))

    ρss = steadystate(H, c_ops)
    corr = correlation_2op_1t(H, ρss, τ_l, A, B, c_ops; kwargs...)

    S = fftshift(fft(corr)) / length(τ_l)

    return ω_l, 2 .* real.(S)
end

function _spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_max::Real, Nsamples::Integer,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector,
    solver::ExponentialSeries;
    kwargs...) where {T1,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    
    (H.dims == A.dims == B.dims) || throw(DimensionMismatch("The dimensions of H, A and B must be the same"))
    Hdims = prod(H.dims)

    L = liouvillian(H, c_ops)

    ω_max = abs(ω_max)
    dω = 2*ω_max/(Nsamples-1)
    ω_l = -ω_max:dω:ω_max

    vals, vecs = eigen(L)
    # The steadystate should be always the last eigenvector
    ρss = QuantumObject(reshape(vecs[:,end], Hdims...), dims=A.dims)
    ρss /= tr(ρss)
    ρss = (ρss + ρss') / 2 # Make sure it's hermitian
    ρ0 = B.data * ρss.data
    v = vecs \ reshape(ρ0, :)
    
    amps = map(i->tr(A.data * reshape(v[i] * vecs[:,i], Hdims...)), 1:length(vals))
    rates = vals
    push!(amps, -expect(A, ρss) * expect(B, ρss))
    push!(rates, 0)
    idxs = sortperm(rates, by=abs)
    # For stability reasons, we take the modulus of the amplitudes
    amps, rates = abs.(amps[idxs]), rates[idxs]
    idxs = amps .> solver.tol
    amps, rates = amps[idxs], rates[idxs]
    
    spec = real.(transpose(amps) * (2 ./ (transpose(1im .* ω_l) .- rates)))[1,:]
    # spec2 = map(ω->real(sum(amps .* (2 ./ (1im * ω .- rates)))), ω_l)

    return ω_l, spec
end