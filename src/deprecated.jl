#=
This file gathers all the deprecated names (structures, functions, or variables) which will be removed in the future major release.

- Before the major release, the deprecated names will just throw errors when they are called.
- If the deprecated names were once exported, we will still export them here until next major release.
- If we decide to push a major release, cleanup this file.

Example:

export deprecated_foo

function deprecated_foo(args...; kwargs...)
    error("`deprecated_foo` has been deprecated and will be removed in next major release, please use `new_foo` instead.")
end
=#

export FFTCorrelation

FFTCorrelation() = error(
    "`FFTCorrelation` has been deprecated and will be removed in next major release, please use `spectrum_correlation_fft` to calculate the spectrum with FFT method instead.",
)

correlation_3op_2t(
    H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector,
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    C::QuantumObject{<:AbstractArray{T5},OperatorQuantumObject},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    kwargs...,
) where {
    T1,
    T2,
    T3,
    T4,
    T5,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
} = error(
    "The parameter order of `correlation_3op_2t` has been changed, please use `?correlation_3op_2t` to check the updated docstring.",
)

correlation_3op_1t(
    H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    C::QuantumObject{<:AbstractArray{T5},OperatorQuantumObject},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    kwargs...,
) where {
    T1,
    T2,
    T3,
    T4,
    T5,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
} = error(
    "The parameter order of `correlation_3op_1t` has been changed, please use `?correlation_3op_1t` to check the updated docstring.",
)

correlation_2op_2t(
    H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    t_l::AbstractVector,
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    reverse::Bool = false,
    kwargs...,
) where {
    T1,
    T2,
    T3,
    T4,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
} = error(
    "The parameter order of `correlation_2op_2t` has been changed, please use `?correlation_2op_2t` to check the updated docstring.",
)

correlation_2op_1t(
    H::QuantumObject{<:AbstractArray{T1},HOpType},
    ψ0::QuantumObject{<:AbstractArray{T2},StateOpType},
    τ_l::AbstractVector,
    A::QuantumObject{<:AbstractArray{T3},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T4},OperatorQuantumObject},
    c_ops::Union{Nothing,AbstractVector,Tuple} = nothing;
    reverse::Bool = false,
    kwargs...,
) where {
    T1,
    T2,
    T3,
    T4,
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
} = error(
    "The parameter order of `correlation_2op_1t` has been changed, please use `?correlation_2op_1t` to check the updated docstring.",
)
