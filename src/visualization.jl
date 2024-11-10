export plot_wigner

plot_wigner(
    state::QuantumObject{<:AbstractArray{T},OpType};
    library::Union{Val,Symbol} = Val(:CairoMakie),
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}} =
    plot_wigner(makeVal(library), state; kwargs...)
