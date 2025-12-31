@doc raw"""
    matrix_histogram(
        library::Val{:Makie},
        M::Union{QuantumObject,AbstractMatrix};
        method::Union{Symbol,Val} = Val(:real),
        xbasis::Union{Nothing,AbstractVector}=nothing,
        ybasis::Union{Nothing,AbstractVector}=nothing,
        limits::Union{Nothing,Tuple{Real,Real}}=nothing,
        azimuth::Real = -60,
        elevation::Real = 45,
        bars_spacing::Real = 0.2,
        shading::Bool = true,
        colorbar::Bool = true,
        colorbar_label::Union{Nothing,AbstractString} = nothing,
        colormap = :viridis,
        colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
        location::Union{GridPosition,Nothing} = nothing,
        kwargs...,
    )

Plot a 3D histogram for the elements of matrix `M`.

# Arguments
- `library::Val{:Makie}`: The plotting library to use.
- `M::Union{QuantumObject,AbstractMatrix}`: The [`QuantumObject`](@ref) or `AbstractMatrix` for which to be plotted. If it is a [`QuantumObject`](@ref), tt can be either a [`Operator`](@ref) or [`SuperOperator`](@ref).
- `method::Union{Symbol,Val}`: Method to use for plotting the matrix elements. Can be either `:real`, `:imag`, `:abs`, or `:angle`. Default is `Val(:real)`.
- `xbasis::Union{Nothing,AbstractVector}`: List of `x`-ticklabels. If `nothing`, the labels of computational basis will be automatically created. Default is `nothing`.
- `ybasis::Union{Nothing,AbstractVector}`: List of `y`-ticklabels. If `nothing`, the labels of computational basis will be automatically created. Default is `nothing`.
- `limits::Union{Nothing,Tuple{Real,Real}}`: The `z`-axis limits `(z_min, z_max)`. If `nothing`, the limits will be automatically decided. Default is `nothing`.
- `azimuth::Real`: The azimuthal viewing angle (in degrees). Default is `-60`.
- `elevation::Real`: The elevation viewing angle (in degrees). Default is `45`.
- `bars_spacing::Real`: Spacing between bars. Default is `0.2`.
- `shading::Bool`: Whether to shade the dark sides of the bars. Note that the shading is relative to plot's source of light. Default to `true`.
- `colorbar::Bool`: Whether to show the color bar. Default is `true`.
- `colorbar_label::Union{Nothing,AbstractString}`: Label of color bar. If nothing, it will automatically use the string of `method`. Default is `nothing`.
- `colormap`: The color map to use when plotting. Default is `:viridis`.
- `colorrange::Union{Nothing,Tuple{Real,Real}}`: The color bar limits `(min, max)`. If nothing, it will automatically use (`z`-axis) `limits`. Default is `nothing`.
- `location::Union{GridPosition,Nothing}`: The location of the plot in the layout. If `nothing`, the plot is created in a new figure. Default is `nothing`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function.

# Returns
- `fig`: The figure object.
- `ax`: The axis object.
- `ms`: The meshscatter object.

!!! note "Import library first"
    [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) must first be imported before using this function. This can be done by importing one of the available backends, such as [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie), [`GLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/GLMakie), or [`WGLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/WGLMakie).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `method = Val(:real)`, `Val(:imag)`, `Val(:abs)`, and `Val(:angle)` instead of `:real`, `:imag`, `:abs`, and `:angle`, respectively. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function QuantumToolbox.matrix_histogram(
    library::Val{:Makie},
    M::Union{QuantumObject{QT},AbstractMatrix{MT}};
    method::Union{Symbol,Val} = Val(:real),
    xbasis::Union{Nothing,AbstractVector} = nothing,
    ybasis::Union{Nothing,AbstractVector} = nothing,
    limits::Union{Nothing,Tuple{Real,Real}} = nothing,
    azimuth::Real = -60,
    elevation::Real = 45,
    bars_spacing::Real = 0.2,
    shading::Bool = true,
    colorbar::Bool = true,
    colorbar_label::Union{Nothing,AbstractString} = nothing,
    colormap = :viridis,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    location::Union{GridPosition,Nothing} = nothing,
    kwargs...,
) where {QT<:Union{Operator,SuperOperator},MT<:Number}
    fig, location = _getFigAndLocation(location)
    lyt = GridLayout(location)

    # generate x and y data
    # y-axis means row, x-axis means column
    Ny, Nx = size(M)
    xdata = 0:(Nx-1)   # computational basis
    ydata = 0:(Ny-1)   # computational basis

    # generate z data
    ## z0 mean the origin of all bars
    ## the matrix data (zdata) will be handled by scaling each marker size in meshscatter! later
    z0 = zeros(Nx, Ny)
    zdata = vec(_handle_matrix_plot_data(M, makeVal(method)))

    # handle x and y ticks
    xbasis = isnothing(xbasis) ? _gen_default_bra_labels(M, xdata) : xbasis
    ybasis = isnothing(ybasis) ? _gen_default_ket_labels(M, ydata) : ybasis
    length(xbasis) == Nx ||
        throw(ArgumentError("Length of xbasis ($(length(xbasis))) does not match matrix size: ($Nx)"))
    length(ybasis) == Ny ||
        throw(ArgumentError("Length of ybasis ($(length(ybasis))) does not match matrix size: ($Ny)"))
    xticks = (xdata, xbasis)
    yticks = (ydata, ybasis)

    # handle z-axis limits and colorrange
    if isnothing(limits)
        z_min = minimum(zdata)
        z_max = maximum(zdata)

        # must include 0 between the interval: [z_min, z_max]
        z_min = z_min > 0 ? 0 : z_min
        z_max = z_max < 0 ? 0 : z_max
        limits = (z_min, z_max)
    end
    colorrange = isnothing(colorrange) ? limits : colorrange

    # handle bar's edge length in x and y directions
    (0 ≤ bars_spacing ≤ 1) ||
        throw(ArgumentError("Invalid keyword argument bars_spacing = $(bars_spacing), should be: 0 ≤ bars_spacing ≤ 1"))
    bar_edge_length = 1 - bars_spacing

    # create 3D axis
    ax = Axis3(
        lyt[1, 1];
        aspect = (Nx / Ny, 1, 0.5), # this ratio makes the elements always be squared size (even in non-squared matrix cases)
        limits = ((-0.5, Nx - 0.5), (-0.5, Ny - 0.5), limits), # -0.5 in x and y originates from the half width of bar
        azimuth = deg2rad(azimuth),
        elevation = deg2rad(elevation),
        xticks = xticks,
        yticks = yticks,
        yreversed = true,
        xlabelvisible = false,
        ylabelvisible = false,
        zlabelvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xgridvisible = false,
        ygridvisible = false,
    )

    # plot meshscatter
    ms = meshscatter!(
        ax,
        xdata,
        ydata,
        z0;
        marker = Rect3f((-0.5, -0.5, 0), (1, 1, 1)), # corner position of the bar (x=-0.5,y=-0.5,z=0); original width (1,1,1)
        markersize = Vec3f.(bar_edge_length, bar_edge_length, zdata), # scale the size of the marker
        shading = shading,
        color = zdata,
        colormap = colormap,
        colorrange = colorrange,
        kwargs...,
    )

    # handle color bar
    if colorbar
        label = isnothing(colorbar_label) ? string(getVal(method)) : colorbar_label
        Colorbar(lyt[1, 2], ms, label = label)
    end

    return fig, ax, ms
end
