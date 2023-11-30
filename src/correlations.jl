abstract type SpectrumSolver end
struct FFTCorrelation <: SpectrumSolver end
struct ExponentialSeries <: SpectrumSolver
    tol::Real
    calc_steadystate::Bool
end

ExponentialSeries(;tol=1e-14,calc_steadystate=false) = ExponentialSeries(tol,calc_steadystate)

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

    kwargs2 = (;kwargs...)
    kwargs2 = merge(kwargs2, (saveat = collect(t_l),))
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
    ω_list::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    solver::MySolver=ExponentialSeries(),
        kwargs...)

Returns the emission spectrum ``S(\omega) = \int_{-\infty}^\infty \expval{\hat{A}(\tau) \hat{B}(0)} e^{-i \omega \tau} d \tau``.
"""
function spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_list::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    solver::MySolver=ExponentialSeries(),
    kwargs...) where {T1,T2,T3,
            HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}, MySolver<:SpectrumSolver}
    
    return _spectrum(H, ω_list, A, B, c_ops, solver; kwargs...)
end

function _spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_list::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector,
    solver::FFTCorrelation;
    kwargs...) where {T1,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    
    Nsamples = length(ω_list)
    ω_max = abs(maximum(ω_list))
    dω = 2*ω_max/(Nsamples-1)
    ω_l = -ω_max:dω:ω_max

    T = 2π/(ω_l[2]-ω_l[1])
    τ_l = range(0, T, length=length(ω_l))

    ρss = steadystate(H, c_ops)
    corr = correlation_2op_1t(H, ρss, τ_l, A, B, c_ops; kwargs...)

    S = fftshift(fft(corr)) / length(τ_l)

    return ω_l, 2 .* real.(S)
end

function _spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_list::AbstractVector,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector,
    solver::ExponentialSeries;
    kwargs...) where {T1,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    
    (H.dims == A.dims == B.dims) || throw(DimensionMismatch("The dimensions of H, A and B must be the same"))

    L = liouvillian(H, c_ops)

    ω_l = ω_list

    rates, vecs = eigen(L)

    # Get the steady state and update the corresponding vector
    # ss_idx = argmin(abs2.(rates))
    # ρss = vec2mat(@view(vecs[:,ss_idx]))
    # ρss2 = (ρss + ρss') / 2
    # ρss2 ./= tr(ρss2)
    # ρss .= ρss2
    # ρss = steadystate(L).data # This solves the problem when multiple states have zero eigenvalue
    if solver.calc_steadystate
        ρss = steadystate(L).data
    else
        ss_idx = argmin(abs2.(rates))
        ρss = vec2mat(@view(vecs[:,ss_idx]))
        ρss2 = (ρss + ρss') / 2
        ρss2 ./= tr(ρss2)
        ρss .= ρss2
    end
    
    ρ0 = B.data * ρss
    v = vecs \ mat2vec(ρ0)
    
    amps = map(i->v[i] * tr(A.data * vec2mat(@view(vecs[:,i]))), eachindex(rates))
    idxs = findall(x -> abs(x) > solver.tol, amps)
    amps, rates = amps[idxs], rates[idxs]
    # @. amps = abs(amps)
    # idxs = findall(x -> real(x) < 0, amps)
    # @. amps[idxs] -= 2*real(amps[idxs])
    
    spec = map(ω->2*real(sum(@. amps * (1 / (1im * ω - rates)))), ω_l)

    return ω_l, spec
end
