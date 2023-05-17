function _correlation_2t(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector, τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    C::QuantumObject{<:AbstractArray{T5},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    alg=Vern7(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks=[],
    kwargs...) where {T1,T2,T3,T4,T5,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    (H.dims == ψ0.dims && H.dims == A.dims &&
    H.dims == B.dims && H.dims == C.dims) || throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    kwargs2 = merge(kwargs, Dict(:saveat => collect(t_l)))
    ρt = mesolve(H, ψ0, t_l, c_ops; alg=alg, H_t=H_t, params=params, progress=progress, callbacks=callbacks, kwargs2...).states

    corr = map((t,ρ)->mesolve(H, C*ρ*A, τ_l .+ t, c_ops; e_ops=[B], alg=alg, H_t=H_t, params=params,
                        progress=progress, callbacks=callbacks, kwargs...).expect[1,:], t_l, ρt)

    corr
end

@doc raw"""
    correlation_2op_1t(H::QuantumObject,
        ψ0::QuantumObject,
        τ_l::AbstractVector,
        A::QuantumObject,
        B::QuantumObject,
        c_ops::AbstractVector=[];
        alg=Vern7(),
        H_t=nothing,
        params::AbstractVector=[],
        progress::Bool=true,
        callbacks=[],
        reverse::Bool=false,
        kwargs...)

Returns the correlation function of two operators `\hat{A}` and `\hat{B}` at different times ``\expval{\hat{A}(\tau) \hat{B}(0)}``.
"""
function correlation_2op_1t(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    alg=Vern7(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks=[],
    reverse::Bool=false,
    kwargs...) where {T1,T2,T3,T4,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject}}
    
    (H.dims == ψ0.dims && H.dims == A.dims &&
    H.dims == B.dims) || throw(ErrorException("The two operators are not of the same Hilbert dimension."))

    C = eye(prod(H.dims), dims=H.dims)
    if reverse
        corr = _correlation_2t(H, ψ0, [0], τ_l, A, B, C, c_ops; alg=alg, H_t=H_t, 
        params=params, progress=progress, callbacks=callbacks, kwargs...)
    else
        corr = _correlation_2t(H, ψ0, [0], τ_l, C, A, B, c_ops; alg=alg, H_t=H_t, 
        params=params, progress=progress, callbacks=callbacks, kwargs...)
    end
    
    corr[1]
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

Returns the emission spectrum `S(\omega) = \int_{-\infty}^\infty \expval{\hat{A}(\tau) \hat{B}(0)} e^{-i \omega \tau} d \tau`.
"""
function spectrum(H::QuantumObject{<:AbstractArray{T1},HOpType},
    ω_max::Real, Nsamples::Integer,
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    c_ops::AbstractVector=[];
    alg=Vern7(),
    H_t=nothing,
    params::AbstractVector=[],
    progress::Bool=true,
    callbacks=[],
    kwargs...) where {T1,T2,T3,HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    
    ω_max = abs(ω_max)
    dω = 2*ω_max/(Nsamples-1)
    ω_l = -ω_max:dω:ω_max

    T = 2π/(ω_l[2]-ω_l[1])
    τ_l = range(0,T, length=length(ω_l))

    ρss = steadystate(H, c_ops)
    corr = correlation_2op_1t(H, ρss, τ_l, A, B, c_ops; alg=alg, H_t=H_t,
    params=params, progress=progress, callbacks=callbacks, kwargs...)

    S = fftshift(fft(corr)) / length(τ_l)

    return ω_l, 2 .* real.(S)
end