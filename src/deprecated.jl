#=
This file gathers all the deprecated names (structures, functions, or variables) which will be removed in the future major release.

- Before the major release, the deprecated names will show warnings or just throw errors when they are called.
- If the deprecated names were once exported, we will still export them here until next major release.
- If we decide to push a major release, cleanup this file.

Example 1 [throw errors if the deprecation is fundamental, and there is no replacement for it]:
```
export deprecated_foo
deprecated_foo(args...; kwargs...) = error("`deprecated_foo` is deprecated and will be removed in next major release.")
```

Example 2 ["force" show warning and tell the users that there is a replacement for the deprecated function]:

```
export deprecated_foo
function deprecated_foo(args...; kwargs...)
    Base.depwarn("`deprecated_foo` is deprecated and will be removed in next major release, use `new_foo` instead.", :deprecated_foo, force = true)
    new_foo(args...; kwargs...)
end
```
=#

export FFTCorrelation
FFTCorrelation() = error(
    "`FFTCorrelation` is deprecated and will be removed in next major release, use `spectrum_correlation_fft` to calculate the spectrum with FFT method instead.",
)

function correlation_3op_2t(
        H::QuantumObject{HOpType},
        ψ0::QuantumObject{StateOpType},
        t_l::AbstractVector,
        τ_l::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    Base.depwarn(
        "The argument order of `correlation_3op_2t(H, ψ0, t_l, τ_l, A, B, C, c_ops)` is changed to `correlation_3op_2t(H, ψ0, t_l, τ_l, c_ops, A, B, C)`, use `?correlation_3op_2t` to check the updated docstring.",
        :correlation_3op_2t,
        force = true,
    )
    return correlation_3op_2t(H, ψ0, t_l, τ_l, c_ops, A, B, C; kwargs...)
end

function correlation_3op_1t(
        H::QuantumObject{HOpType},
        ψ0::QuantumObject{StateOpType},
        τ_l::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    Base.depwarn(
        "The argument order of `correlation_3op_1t(H, ψ0, τ_l, A, B, C, c_ops)` is changed to `correlation_3op_1t(H, ψ0, τ_l, c_ops, A, B, C)`, use `?correlation_3op_1t` to check the updated docstring.",
        :correlation_3op_1t,
        force = true,
    )
    return correlation_3op_1t(H, ψ0, τ_l, c_ops, A, B, C; kwargs...)
end

function correlation_2op_2t(
        H::QuantumObject{HOpType},
        ψ0::QuantumObject{StateOpType},
        t_l::AbstractVector,
        τ_l::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        reverse::Bool = false,
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    Base.depwarn(
        "The argument order of `correlation_2op_2t(H, ψ0, t_l, τ_l, A, B, c_ops)` is changed to `correlation_2op_2t(H, ψ0, t_l, τ_l, c_ops, A, B)`, use `?correlation_2op_2t` to check the updated docstring.",
        :correlation_2op_2t,
        force = true,
    )
    return correlation_2op_2t(H, ψ0, t_l, τ_l, c_ops, A, B; reverse = reverse, kwargs...)
end

function correlation_2op_1t(
        H::QuantumObject{HOpType},
        ψ0::QuantumObject{StateOpType},
        τ_l::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        c_ops::Union{Nothing, AbstractVector, Tuple} = nothing;
        reverse::Bool = false,
        kwargs...,
    ) where {HOpType <: Union{Operator, SuperOperator}, StateOpType <: Union{Ket, Operator}}
    Base.depwarn(
        "The argument order of `correlation_2op_1t(H, ψ0, τ_l, A, B, c_ops)` is changed to `correlation_2op_1t(H, ψ0, τ_l, c_ops, A, B)`, use `?correlation_2op_1t` to check the updated docstring.",
        :correlation_2op_1t,
        force = true,
    )
    return correlation_2op_1t(H, ψ0, τ_l, c_ops, A, B; reverse = reverse, kwargs...)
end

export ProgressBar
function ProgressBar(args...; kwargs...)
    # Use error instead of depwarn, since ProgressBar and Progress have different arguments
    return error(
        "`ProgressBar` is deprecated and will be removed in next major release. Use `Progress` from `ProgressMeter.jl` instead.",
    )
end

export liouvillian_generalized
function liouvillian_generalized(args...; kwargs...)
    Base.depwarn(
        "`liouvillian_generalized` is deprecated and will be removed in next major release, use `liouvillian_dressed_nonsecular` instead.",
        :liouvillian_generalized,
        force = true,
    )
    return liouvillian_dressed_nonsecular(args...; kwargs...)
end

export steadystate_floquet
function steadystate_floquet(args...; kwargs...)
    Base.depwarn(
        "`steadystate_floquet` is deprecated and will be removed in next major release, use `steadystate_fourier` instead.",
        :steadystate_floquet,
        force = true,
    )
    return steadystate_fourier(args...; kwargs...)
end
