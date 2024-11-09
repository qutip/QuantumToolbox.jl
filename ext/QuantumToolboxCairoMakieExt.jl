module QuantumToolboxCairoMakieExt

using QuantumToolbox
using CairoMakie

function QuantumToolbox.plot_wigner(
    library::Val{:CairoMakie},
    state::QuantumObject{<:AbstractArray{T},OpType},
    xvec::Union{Nothing,AbstractVector} = nothing,
    yvec::Union{Nothing,AbstractVector} = nothing;
    g::Real = √2,
    method::WignerSolver = WignerClenshaw(),
    projection::String = "2d",
    fig::Union{Figure,Nothing} = nothing,
    ax::Union{Axis,Nothing} = nothing,
    colorbar::Bool = false,
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    projection == "2d" || projection == "3d" || throw(ArgumentError("Unsupported projection: $projection"))

    return _plot_wigner(
        library,
        state,
        xvec,
        yvec,
        Val(Symbol(projection)),
        g,
        method,
        fig,
        ax,
        colorbar;
        kwargs...
    )
end

function _plot_wigner(
    ::Val{:CairoMakie},
    state::QuantumObject{<:AbstractArray{T},OpType},
    xvec::AbstractVector,
    yvec::AbstractVector,
    projection::Val{Symbol("2d")},
    g::Real,
    method::WignerSolver,
    fig::Union{Figure,Nothing},
    ax::Union{Axis,Nothing},
    colorbar::Bool;
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    fig, ax = _getFigAx(fig, ax)

    gridPos = _gridPosFromAx(ax)
    CairoMakie.delete!(ax)

    lyt = GridLayout(gridPos)
    ax = Axis(lyt[1, 1])

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => :RdBu, :colorrange => (-wlim, wlim)), kwargs)
    hm = heatmap!(ax, xvec, yvec, wig; kwargs...)

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
    projection::Val{Symbol("3d")},
    g::Real,
    method::WignerSolver,
    fig::Union{Figure,Nothing},
    ax::Union{Axis,Nothing},
    colorbar::Bool;
    kwargs...,
) where {T,OpType<:Union{BraQuantumObject,KetQuantumObject,OperatorQuantumObject}}
    fig, ax = _getFigAx(fig, ax)

    gridPos = _gridPosFromAx(ax)
    CairoMakie.delete!(ax)

    lyt = GridLayout(gridPos)
    ax = Axis3(lyt[1, 1], azimuth = 1.775pi, elevation = pi / 16, protrusions = (30, 90, 30, 30), viewmode = :stretch)

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => :RdBu, :colorrange => (-wlim, wlim)), kwargs)
    surf = surface!(ax, xvec, yvec, wig; kwargs...)

    if colorbar
        Colorbar(lyt[1, 2], surf)
    end

    ax.xlabel = L"\Re(\alpha)"
    ax.ylabel = L"\Im(\alpha)"
    ax.zlabel = "Wigner function"
    return fig, ax, surf
end

_getFigAx(fig::Figure, ax::Axis) = fig, ax
_getFigAx(fig::Figure, ::Nothing) = fig, Axis(fig[1, 1])
_getFigAx(::Nothing, ax::Axis) = _figFromChildren(ax), ax
function _getFigAx(::Nothing, ::Nothing)
    fig = Figure()
    ax = Axis(fig[1, 1])
    return fig, ax
end

_figFromChildren(children) = _figFromChildren(children.parent)
_figFromChildren(fig::Figure) = fig

function _gridPosFromAx(ax::Axis)
    content = CairoMakie.Makie.GridLayoutBase.gridcontent(ax)
    gl, sp, si = content.parent, content.span, content.side
    return GridPosition(gl, sp, si)
end

end