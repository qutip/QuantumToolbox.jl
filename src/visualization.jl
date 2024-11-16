export plot_wigner

@doc raw"""
    plot_wigner(
        state::QuantumObject{<:AbstractArray{T},OpType}; 
        library::Union{Val,Symbol}=Val(:CairoMakie), 
        kwargs...
    )

Plot the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) of `state` using the [`wigner`](@ref) function. 
    
The `library` keyword argument specifies the plotting library to use, defaulting to `CairoMakie`. Note that plotting libraries must first be imported before using them with this function.

# Arguments
- `state::QuantumObject{<:AbstractArray{T},OpType}`: The quantum state for which to plot the Wigner distribution.
- `library::Union{Val,Symbol}`: The plotting library to use. Default is `Val(:CairoMakie)`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. See the documentation for the specific plotting library for more information.
"""
plot_wigner(
    state::QuantumObject{<:AbstractArray{T},OpType};
    library::Union{Val,Symbol} = Val(:CairoMakie),
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}} =
    plot_wigner(makeVal(library), state; kwargs...)
