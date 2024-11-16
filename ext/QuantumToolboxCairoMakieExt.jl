module QuantumToolboxCairoMakieExt

using QuantumToolbox
using CairoMakie

@doc raw"""
    plot_wigner(
        library::Val{:CairoMakie},
        state::QuantumObject{<:AbstractArray{T},OpType};
        xvec::Union{Nothing,AbstractVector} = nothing,
        yvec::Union{Nothing,AbstractVector} = nothing,
        g::Real = √2,
        method::WignerSolver = WignerClenshaw(),
        projection::Union{Val,Symbol} = Val(:two_dim),
        location::Union{GridPosition,Nothing} = nothing,
        colorbar::Bool = false,
        kwargs...
    )

Plot the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) of `state` using the [CairoMakie](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie) plotting library.

Note that CairoMakie must first be imported before using this function.

# Arguments
- `library::Val{:CairoMakie}`: The plotting library to use.
- `state::QuantumObject`: The quantum state for which the Wigner function is calculated. It can be either a [`KetQuantumObject`](@ref), [`BraQuantumObject`](@ref), or [`OperatorQuantumObject`](@ref).
- `xvec::AbstractVector`: The x-coordinates of the phase space grid. Default is # TODO
- `yvec::AbstractVector`: The y-coordinates of the phase space grid. Default is # TODO
- `g::Real`: The scaling factor related to the value of ``\hbar`` in the commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2``.
- `method::WignerSolver`: The method used to calculate the Wigner function. It can be either `WignerLaguerre()` or `WignerClenshaw()`, with `WignerClenshaw()` as default. The `WignerLaguerre` method has the optional `parallel` and `tol` parameters, with default values `true` and `1e-14`, respectively.
- `projection::Union{Val,Symbol}`: Wheather to plot the Wigner function in 2D or 3D. It can be either `Val(:two_dim)` or `Val(:three_dim)`, with `Val(:two_dim)` as default.
- `location::Union{GridPosition,Nothing}`: The location of the plot in the layout. If `nothing`, the plot is created in a new figure. Default is `nothing`.
- `colorbar::Bool`: Whether to include a colorbar in the plot. Default is `false`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. 

# Returns
- `fig`: The figure object.
- `ax`: The axis object.
- `hm`: Either the heatmap or surface object, depending on the projection.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:two_dim)` and `Val(:three_dim)` instead of `:two_dim` and `:three_dim`, respectively. Also, specify the library as `Val(:CairoMakie)` See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function QuantumToolbox.plot_wigner(
    library::Val{:CairoMakie},
    state::QuantumObject{<:AbstractArray{T},OpType};
    xvec::Union{Nothing,AbstractVector} = nothing,
    yvec::Union{Nothing,AbstractVector} = nothing,
    g::Real = √2,
    method::WignerSolver = WignerClenshaw(),
    projection::Union{Val,Symbol} = Val(:two_dim),
    location::Union{GridPosition,Nothing} = nothing,
    colorbar::Bool = false,
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    QuantumToolbox.getVal(projection) == :two_dim ||
        QuantumToolbox.getVal(projection) == :three_dim ||
        throw(ArgumentError("Unsupported projection: $projection"))

    return _plot_wigner(
        library,
        state,
        xvec,
        yvec,
        QuantumToolbox.makeVal(projection),
        g,
        method,
        location,
        colorbar;
        kwargs...,
    )
end

function _plot_wigner(
    ::Val{:CairoMakie},
    state::QuantumObject{<:AbstractArray{T},OpType},
    xvec::AbstractVector,
    yvec::AbstractVector,
    projection::Val{:two_dim},
    g::Real,
    method::WignerSolver,
    location::Union{GridPosition,Nothing},
    colorbar::Bool;
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    fig, location = _getFigAndLocation(location)

    lyt = GridLayout(location)

    ax = Axis(lyt[1, 1])

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => Reverse(:RdBu), :colorrange => (-wlim, wlim)), kwargs)
    hm = heatmap!(ax, xvec, yvec, wig'; kwargs...)

    if colorbar
        Colorbar(lyt[1, 2], hm)
    end

    ax.xlabel = L"\Re(\alpha)"
    ax.ylabel = L"\Im(\alpha)"
    return fig, ax, hm
end

function _plot_wigner(
    ::Val{:CairoMakie},
    state::QuantumObject{<:AbstractArray{T},OpType},
    xvec::AbstractVector,
    yvec::AbstractVector,
    projection::Val{:three_dim},
    g::Real,
    method::WignerSolver,
    location::Union{GridPosition,Nothing},
    colorbar::Bool;
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    fig, location = _getFigAndLocation(location)

    lyt = GridLayout(location)

    ax = Axis3(lyt[1, 1], azimuth = 1.775pi, elevation = pi / 16, protrusions = (30, 90, 30, 30), viewmode = :stretch)

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => :RdBu, :colorrange => (-wlim, wlim)), kwargs)
    surf = surface!(ax, xvec, yvec, wig'; kwargs...)

    if colorbar
        Colorbar(lyt[1, 2], surf)
    end

    ax.xlabel = L"\Re(\alpha)"
    ax.ylabel = L"\Im(\alpha)"
    ax.zlabel = "Wigner function"
    return fig, ax, surf
end

function _getFigAndLocation(location::Nothing)
    fig = Figure()
    return fig, fig[1, 1]
end
function _getFigAndLocation(location::GridPosition)
    fig = _figFromChildren(location.layout)
    return fig, location
end

_figFromChildren(children) = _figFromChildren(children.parent)
_figFromChildren(fig::Figure) = fig
_figFromChildren(::Nothing) = throw(ArgumentError("No Figure has been found at the top of the layout hierarchy."))

end