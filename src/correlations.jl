export correlation_3op_2t, correlation_3op_1t, correlation_2op_2t, correlation_2op_1t

function _check_correlation_time_list(tlist::AbstractVector)
    any(t -> t == 0, tlist) ||
        throw(ArgumentError("The time list for calculating correlation function must contain the element `0`"))
    all(>=(0), tlist) ||
        throw(ArgumentError("All the elements in the time list for calculating correlation function must be positive."))
    return nothing
end

@doc raw"""
    correlation_3op_2t(H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject,
        C::QuantumObject;
        kwargs...)

Returns the two-times correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

If the initial state `ψ0` is given as `nothing`, then the [`steadystate`](@ref) will be used as the initial state. Note that this is only implemented if `H` is constant ([`QuantumObject`](@ref)).
"""
function correlation_3op_2t(
    H::AbstractQuantumObject{HOpType},
    ψ0::Union{Nothing,QuantumObject{StateOpType}},
    tlist::AbstractVector,
    τlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{OperatorQuantumObject},
    B::QuantumObject{OperatorQuantumObject},
    C::QuantumObject{OperatorQuantumObject};
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    # check tlist and τlist
    _check_correlation_time_list(tlist)
    _check_correlation_time_list(τlist)

    L = liouvillian(H, c_ops)
    if ψ0 isa Nothing
        ψ0 = steadystate(L)
    end

    check_dimensions(L, ψ0, A, B, C)

    kwargs2 = merge((saveat = collect(tlist),), (; kwargs...))
    ρt_list = mesolve(L, ψ0, tlist; kwargs2...).states

    corr = map((t, ρt) -> mesolve(L, C * ρt * A, τlist .+ t, e_ops = [B]; kwargs...).expect[1, :], tlist, ρt_list)

    # make the output correlation Matrix align with QuTiP
    # 1st dimension corresponds to tlist
    # 2nd dimension corresponds to τlist
    return reduce(vcat, transpose.(corr))
end

@doc raw"""
    correlation_3op_1t(H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject,
        C::QuantumObject;
        kwargs...)

Returns the one-time correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

If the initial state `ψ0` is given as `nothing`, then the [`steadystate`](@ref) will be used as the initial state. Note that this is only implemented if `H` is constant ([`QuantumObject`](@ref)).
"""
function correlation_3op_1t(
    H::AbstractQuantumObject{HOpType},
    ψ0::Union{Nothing,QuantumObject{StateOpType}},
    τlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{OperatorQuantumObject},
    B::QuantumObject{OperatorQuantumObject},
    C::QuantumObject{OperatorQuantumObject};
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    corr = correlation_3op_2t(H, ψ0, [0], τlist, c_ops, A, B, C; kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end

@doc raw"""
    correlation_2op_2t(H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        tlist::AbstractVector,
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject;
        reverse::Bool=false,
        kwargs...)

Returns the two-times correlation function of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

If the initial state `ψ0` is given as `nothing`, then the [`steadystate`](@ref) will be used as the initial state. Note that this is only implemented if `H` is constant ([`QuantumObject`](@ref)).

When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
"""
function correlation_2op_2t(
    H::AbstractQuantumObject{HOpType},
    ψ0::Union{Nothing,QuantumObject{StateOpType}},
    tlist::AbstractVector,
    τlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{OperatorQuantumObject},
    B::QuantumObject{OperatorQuantumObject};
    reverse::Bool = false,
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    C = eye(prod(H.dimensions), dims = H.dimensions)
    if reverse
        corr = correlation_3op_2t(H, ψ0, tlist, τlist, c_ops, A, B, C; kwargs...)
    else
        corr = correlation_3op_2t(H, ψ0, tlist, τlist, c_ops, C, A, B; kwargs...)
    end

    return corr
end

@doc raw"""
    correlation_2op_1t(H::AbstractQuantumObject,
        ψ0::Union{Nothing,QuantumObject},
        τlist::AbstractVector,
        c_ops::Union{Nothing,AbstractVector,Tuple},
        A::QuantumObject,
        B::QuantumObject;
        reverse::Bool=false,
        kwargs...)

Returns the one-time correlation function of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` for a given initial state ``|\psi_0\rangle``.

If the initial state `ψ0` is given as `nothing`, then the [`steadystate`](@ref) will be used as the initial state. Note that this is only implemented if `H` is constant ([`QuantumObject`](@ref)).

When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(0) \hat{B}(\tau) \right\rangle``.
"""
function correlation_2op_1t(
    H::AbstractQuantumObject{HOpType},
    ψ0::Union{Nothing,QuantumObject{StateOpType}},
    τlist::AbstractVector,
    c_ops::Union{Nothing,AbstractVector,Tuple},
    A::QuantumObject{OperatorQuantumObject},
    B::QuantumObject{OperatorQuantumObject};
    reverse::Bool = false,
    kwargs...,
) where {
    HOpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    StateOpType<:Union{KetQuantumObject,OperatorQuantumObject},
}
    corr = correlation_2op_2t(H, ψ0, [0], τlist, c_ops, A, B; reverse = reverse, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end
