export plot_wigner

plot_wigner(library::Val{T}, args...; kwargs...) where {T} =
    throw(ArgumentError("Unsupported visualization library: $(getVal(library))"))
plot_wigner(library::Symbol, args...; kwargs...) = plot_wigner(Val(library), args...; kwargs...)