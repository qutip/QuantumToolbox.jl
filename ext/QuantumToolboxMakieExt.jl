module QuantumToolboxMakieExt

using QuantumToolbox
import QuantumToolbox: _state_to_bloch

import LinearAlgebra: cross, deg2rad, normalize, size
import Makie:
    Axis,
    Axis3,
    LScene,
    Colorbar,
    Figure,
    GridLayout,
    heatmap!,
    surface!,
    barplot!,
    GridPosition,
    @L_str,
    Reverse,
    ylims!,
    RGBAf,
    Sphere,
    lines!,
    scatter!,
    arrows3d!,
    text!,
    Point3f,
    mesh!,
    RGBf,
    Vec3f,
    Point3f,
    NoShading,
    cameracontrols,
    update_cam!,
    cam3d!

@doc raw"""
    plot_wigner(
        library::Val{:Makie},
        state::QuantumObject{OpType};
        xvec::Union{Nothing,AbstractVector} = nothing,
        yvec::Union{Nothing,AbstractVector} = nothing,
        g::Real = √2,
        method::WignerSolver = WignerClenshaw(),
        projection::Union{Val,Symbol} = Val(:two_dim),
        location::Union{GridPosition,Nothing} = nothing,
        colorbar::Bool = false,
        kwargs...
    ) where {OpType}

Plot the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) of `state` using the [`Makie`](https://github.com/MakieOrg/Makie.jl) plotting library.

# Arguments
- `library::Val{:Makie}`: The plotting library to use.
- `state::QuantumObject`: The quantum state for which the Wigner function is calculated. It can be either a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `xvec::AbstractVector`: The x-coordinates of the phase space grid. Defaults to a linear range from -7.5 to 7.5 with 200 points.
- `yvec::AbstractVector`: The y-coordinates of the phase space grid. Defaults to a linear range from -7.5 to 7.5 with 200 points.
- `g::Real`: The scaling factor related to the value of ``\hbar`` in the commutation relation ``[x, y] = i \hbar`` via ``\hbar=2/g^2``.
- `method::WignerSolver`: The method used to calculate the Wigner function. It can be either `WignerLaguerre()` or `WignerClenshaw()`, with `WignerClenshaw()` as default. The `WignerLaguerre` method has the optional `parallel` and `tol` parameters, with default values `true` and `1e-14`, respectively.
- `projection::Union{Val,Symbol}`: Whether to plot the Wigner function in 2D or 3D. It can be either `Val(:two_dim)` or `Val(:three_dim)`, with `Val(:two_dim)` as default.
- `location::Union{GridPosition,Nothing}`: The location of the plot in the layout. If `nothing`, the plot is created in a new figure. Default is `nothing`.
- `colorbar::Bool`: Whether to include a colorbar in the plot. Default is `false`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. 

# Returns
- `fig`: The figure object.
- `ax`: The axis object.
- `hm`: Either the heatmap or surface object, depending on the projection.

!!! note "Import library first"
    [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) must first be imported before using this function. This can be done by importing one of the available backends, such as [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie), [`GLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/GLMakie), or [`WGLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/WGLMakie).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:two_dim)` and `Val(:three_dim)` instead of `:two_dim` and `:three_dim`, respectively. Also, specify the library as `Val(:Makie)` See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function QuantumToolbox.plot_wigner(
    library::Val{:Makie},
    state::QuantumObject{OpType};
    xvec::Union{Nothing,AbstractVector} = LinRange(-7.5, 7.5, 200),
    yvec::Union{Nothing,AbstractVector} = LinRange(-7.5, 7.5, 200),
    g::Real = √2,
    method::WignerSolver = WignerClenshaw(),
    projection::Union{Val,Symbol} = Val(:two_dim),
    location::Union{GridPosition,Nothing} = nothing,
    colorbar::Bool = false,
    kwargs...,
) where {OpType<:Union{Bra,Ket,Operator}}
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
    ::Val{:Makie},
    state::QuantumObject{OpType},
    xvec::AbstractVector,
    yvec::AbstractVector,
    projection::Val{:two_dim},
    g::Real,
    method::WignerSolver,
    location::Union{GridPosition,Nothing},
    colorbar::Bool;
    kwargs...,
) where {OpType<:Union{Bra,Ket,Operator}}
    fig, location = _getFigAndLocation(location)

    lyt = GridLayout(location)

    ax = Axis(lyt[1, 1])

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => Reverse(:RdBu), :colorrange => (-wlim, wlim)), kwargs)
    hm = heatmap!(ax, xvec, yvec, transpose(wig); kwargs...)

    if colorbar
        Colorbar(lyt[1, 2], hm)
    end

    ax.xlabel = L"\textrm{Re}(\alpha)"
    ax.ylabel = L"\textrm{Im}(\alpha)"
    return fig, ax, hm
end

function _plot_wigner(
    ::Val{:Makie},
    state::QuantumObject{OpType},
    xvec::AbstractVector,
    yvec::AbstractVector,
    projection::Val{:three_dim},
    g::Real,
    method::WignerSolver,
    location::Union{GridPosition,Nothing},
    colorbar::Bool;
    kwargs...,
) where {OpType<:Union{Bra,Ket,Operator}}
    fig, location = _getFigAndLocation(location)

    lyt = GridLayout(location)

    ax = Axis3(lyt[1, 1], azimuth = 1.775pi, elevation = pi / 16, protrusions = (30, 90, 30, 30), viewmode = :stretch)

    wig = wigner(state, xvec, yvec; g = g, method = method)
    wlim = maximum(abs, wig)

    kwargs = merge(Dict(:colormap => :RdBu, :colorrange => (-wlim, wlim)), kwargs)
    surf = surface!(ax, xvec, yvec, transpose(wig); kwargs...)

    if colorbar
        Colorbar(lyt[1, 2], surf)
    end

    ax.xlabel = L"\textrm{Re}(\alpha)"
    ax.ylabel = L"\textrm{Im}(\alpha)"
    ax.zlabel = "Wigner function"
    return fig, ax, surf
end

@doc raw"""
    plot_fock_distribution(
        library::Val{:Makie},
        ρ::QuantumObject{SType};
        fock_numbers::Union{Nothing, AbstractVector} = nothing,
        unit_y_range::Bool = true,
        location::Union{GridPosition,Nothing} = nothing,
        kwargs...
    ) where {SType<:Union{Ket,Operator}}

Plot the [Fock state](https://en.wikipedia.org/wiki/Fock_state) distribution of `ρ`. 

# Arguments
- `library::Val{:Makie}`: The plotting library to use.
- `ρ::QuantumObject`: The quantum state for which the Fock state distribution is to be plotted. It can be either a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `location::Union{GridPosition,Nothing}`: The location of the plot in the layout. If `nothing`, the plot is created in a new figure. Default is `nothing`.
- `fock_numbers::Union{Nothing, AbstractVector}`: list of x ticklabels to represent fock numbers, default is `nothing`.
- `unit_y_range::Bool`: Set y-axis limits [0, 1] or not, default is `true`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. 

# Returns
- `fig`: The figure object.
- `ax`: The axis object.
- `bp`: The barplot object.

!!! note "Import library first"
    [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) must first be imported before using this function. This can be done by importing one of the available backends, such as [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie), [`GLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/GLMakie), or [`WGLMakie.jl`](https://github.com/MakieOrg/Makie.jl/tree/master/WGLMakie).

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:two_dim)` and `Val(:three_dim)` instead of `:two_dim` and `:three_dim`, respectively. Also, specify the library as `Val(:Makie)` See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
function QuantumToolbox.plot_fock_distribution(
    library::Val{:Makie},
    ρ::QuantumObject{SType};
    fock_numbers::Union{Nothing,AbstractVector} = nothing,
    unit_y_range::Bool = true,
    location::Union{GridPosition,Nothing} = nothing,
    kwargs...,
) where {SType<:Union{Bra,Ket,Operator}}
    return _plot_fock_distribution(
        library,
        ρ;
        fock_numbers = fock_numbers,
        unit_y_range = unit_y_range,
        location = location,
        kwargs...,
    )
end

function _plot_fock_distribution(
    ::Val{:Makie},
    ρ::QuantumObject{SType};
    fock_numbers::Union{Nothing,AbstractVector} = nothing,
    unit_y_range::Bool = true,
    location::Union{GridPosition,Nothing} = nothing,
    kwargs...,
) where {SType<:Union{Bra,Ket,Operator}}
    ρ = ket2dm(ρ)
    D = prod(ρ.dims)
    isapprox(tr(ρ), 1, atol = 1e-4) || (@warn "The input ρ should be normalized.")

    xvec = 0:(D-1)
    isnothing(fock_numbers) && (fock_numbers = string.(collect(xvec)))

    fig, location = _getFigAndLocation(location)
    lyt = GridLayout(location)
    ax = Axis(lyt[1, 1])

    bp = barplot!(ax, xvec, real(diag(ρ)); kwargs...)

    ax.xticks = (xvec, fock_numbers)
    ax.xlabel = "Fock number"
    ax.ylabel = "Occupation probability"
    unit_y_range && ylims!(ax, 0, 1)

    return fig, ax, bp
end

raw"""
    _getFigAndLocation(location::Nothing)
    
    Create a new figure and return it, together with the GridPosition object pointing to the first cell.

    # Arguments
    - `location::Nothing`

    # Returns
    - `fig`: The figure object.
    - `location`: The GridPosition object pointing to the first cell.
"""
function _getFigAndLocation(location::Nothing)
    fig = Figure()
    return fig, fig[1, 1]
end

raw"""
    _getFigAndLocation(location::GridPosition)
    
    Compute which figure does the location belong to and return it, together with the location itself.

    # Arguments
    - `location::GridPosition`

    # Returns
    - `fig`: The figure object.
    - `location`: The GridPosition object.
"""
function _getFigAndLocation(location::GridPosition)
    fig = _figFromChildren(location.layout)
    return fig, location
end

raw"""
    _figFromChildren(children::GridLayout)

    Recursively find the figure object from the children layout.

    # Arguments
    - `children::GridLayout`

    # Returns
    - Union{Nothing, Figure, GridLayout}: The children's parent object.
"""
_figFromChildren(children) = _figFromChildren(children.parent)

raw"""
    _figFromChildren(fig::Figure)

    Return the figure object

    # Arguments
    - `fig::Figure`

    # Returns
    - `fig`: The figure object.
"""
_figFromChildren(fig::Figure) = fig

raw"""
    _figFromChildren(::Nothing)

    Throw an error if no figure has been found.

    # Arguments
    - `::Nothing`

    # Throws
    - `ArgumentError`: If no figure has been found.
"""
_figFromChildren(::Nothing) = throw(ArgumentError("No Figure has been found at the top of the layout hierarchy."))

function _render_bloch_makie(bloch_vec::Vector{Float64}; location = nothing, kwargs...)
    b = Bloch()
    add_vectors!(b, bloch_vec)
    fig, location = _getFigAndLocation(location)
    fig, ax = render(b; location = location, kwargs...)
    return fig, ax
end

@doc raw"""
    render(b::Bloch; location=nothing)

Render the Bloch sphere visualization from the given [`Bloch`](@ref) object `b`.

# Arguments

- `b::Bloch`: The Bloch sphere object containing states, vectors, and settings to visualize.
- `location::Union{GridPosition,LScene,Nothing}`: The location of the plot in the layout, or `Makie.LScene`. Default is `nothing`.


# Returns

- A tuple `(fig, lscene)` where `fig` is the figure object and `lscene` is the LScene object used for plotting. These can be further manipulated or saved by the user.

# Notes

The keyword argument `location` can be in the either type: 

- `Nothing` (default): Create a new figure and plot the Bloch sphere.
- `GridPosition`: Plot the Bloch sphere in the specified location of the plot in the layout.
- `LScene`: Update the existing Bloch sphere using new data and settings in `b::Bloch` without creating new `Figure` and `LScene` (efficient for drawing animation).
"""
function QuantumToolbox.render(b::Bloch; location = nothing)
    fig, lscene = _setup_bloch_plot!(location)
    _setup_bloch_camara!(b, lscene)
    _draw_bloch_sphere!(b, lscene)
    _add_labels!(b, lscene)

    # plot data fields in Bloch
    _plot_vectors!(b, lscene)
    _plot_lines!(b, lscene)
    _plot_arcs!(b, lscene)
    _plot_points!(b, lscene) # plot points at the end so that they will be on the very top (front) figure layer.

    return fig, lscene
end

raw"""
    _setup_bloch_plot!(location) -> (fig, lscene)

Initialize the Figure and LScene for Bloch sphere visualization.

# Arguments
- `location`: Figure layout position specification, or directly `Makie.LScene` for updating Bloch sphere.

# Returns
- `fig`: Created Makie figure
- `lscene`: Configured LScene object
"""
function _setup_bloch_plot!(location)
    fig, location = _getFigAndLocation(location)
    lscene = LScene(location, show_axis = false, scenekw = (clear = true,))
    return fig, lscene
end

function _setup_bloch_plot!(lscene::LScene)
    # this function only removes all existing Plots in lscene
    # it is useful for users to just update Bloch sphere without creating new figure and lscene (efficient for drawing animation)
    fig = lscene.parent
    empty!(lscene.scene.plots)
    return fig, lscene
end

raw"""
    _setup_bloch_camara!(b::Bloch, lscene)

Setup the distance and view angle of the camara.
"""
function _setup_bloch_camara!(b::Bloch, lscene)
    length(b.view) == 2 || throw(ArgumentError("The length of `Bloch.view` must be 2."))
    cam3d!(lscene.scene, center = false)
    cam = cameracontrols(lscene)
    cam.fov[] = 12 # Set field of view to 12 degrees
    dist = 12      # Set distance from the camera to the Bloch sphere
    update_cam!(lscene.scene, cam, deg2rad(b.view[1]), deg2rad(b.view[2]), dist)
    return nothing
end

raw"""
    _draw_bloch_sphere!(b::Bloch, lscene)

Draw the translucent sphere, axes, and reference circles representing the Bloch sphere surface.
"""
function _draw_bloch_sphere!(b::Bloch, lscene)
    radius = 1.0f0
    sphere_mesh = Sphere(Point3f(0), radius)
    mesh!(
        lscene,
        sphere_mesh;
        color = b.sphere_color,
        alpha = b.sphere_alpha,
        shading = NoShading,
        transparency = true,
        rasterize = 3,
    )

    # X, Y, and Z axes
    axes = [
        [Point3f(1.0, 0, 0), Point3f(-1.0, 0, 0)],  # X-axis
        [Point3f(0, 1.0, 0), Point3f(0, -1.0, 0)],  # Y-axis
        [Point3f(0, 0, 1.0), Point3f(0, 0, -1.0)],  # Z-axis
    ]
    for points in axes
        lines!(lscene, points; color = b.frame_color)
    end

    # highlight circles for XY and XZ planes
    φ = range(0, 2π, length = 100)
    lines!(lscene, [Point3f(cos(φi), sin(φi), 0) for φi in φ]; color = b.frame_color, linewidth = b.frame_width) # XY
    lines!(lscene, [Point3f(cos(φi), 0, sin(φi)) for φi in φ]; color = b.frame_color, linewidth = b.frame_width) # XZ

    # other curves of longitude (with polar angle φ and azimuthal angle θ)
    φ_curve = range(0, 2π, 600)
    θ_vals = [1, 2, 3] * π / 4
    for θi in θ_vals
        x_line = radius * sin.(φ_curve) .* cos(θi)
        y_line = radius * sin.(φ_curve) .* sin(θi)
        z_line = radius * cos.(φ_curve)
        lines!(lscene, x_line, y_line, z_line; color = b.frame_color, alpha = b.frame_alpha, linewidth = b.frame_width)
    end

    # other curves of latitude (with polar angle φ and azimuthal angle θ)
    φ_vals = [1, 3] * π / 4  # missing `2` because XY plane has already be handled above
    θ_curve = range(0, 2π, 600)
    for ϕ in φ_vals
        x_ring = radius * sin(ϕ) .* cos.(θ_curve)
        y_ring = radius * sin(ϕ) .* sin.(θ_curve)
        z_ring = fill(radius * cos(ϕ), length(θ_curve))
        lines!(lscene, x_ring, y_ring, z_ring; color = b.frame_color, alpha = b.frame_alpha, linewidth = b.frame_width)
    end
    return nothing
end

raw"""
    _add_labels!(b::Bloch, lscene)

Add axis labels and state labels to the Bloch sphere.

# Arguments
- `lscene`: LScene object for text placement

Positions standard labels `(x, y, |0⟩, |1⟩)` at appropriate locations.
"""
function _add_labels!(b::Bloch, lscene)
    length(b.xlabel) == 2 || throw(ArgumentError("The length of `Bloch.xlabel` must be 2."))
    length(b.ylabel) == 2 || throw(ArgumentError("The length of `Bloch.ylabel` must be 2."))
    length(b.zlabel) == 2 || throw(ArgumentError("The length of `Bloch.zlabel` must be 2."))
    length(b.xlpos) == 2 || throw(ArgumentError("The length of `Bloch.xlpos` must be 2."))
    length(b.ylpos) == 2 || throw(ArgumentError("The length of `Bloch.ylpos` must be 2."))
    length(b.zlpos) == 2 || throw(ArgumentError("The length of `Bloch.zlpos` must be 2."))

    label_color = parse(RGBf, b.font_color)
    label_size = b.font_size

    (b.xlabel[1] == "") || text!(
        lscene,
        b.xlabel[1],
        position = Point3f(b.xlpos[1], 0, 0),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    (b.xlabel[2] == "") || text!(
        lscene,
        b.xlabel[2],
        position = Point3f(b.xlpos[2], 0, 0),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    (b.ylabel[1] == "") || text!(
        lscene,
        b.ylabel[1],
        position = Point3f(0, b.ylpos[1], 0),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    (b.ylabel[2] == "") || text!(
        lscene,
        b.ylabel[2],
        position = Point3f(0, b.ylpos[2], 0),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    (b.zlabel[1] == "") || text!(
        lscene,
        b.zlabel[1],
        position = Point3f(0, 0, b.zlpos[1]),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    (b.zlabel[2] == "") || text!(
        lscene,
        b.zlabel[2],
        position = Point3f(0, 0, b.zlpos[2]),
        color = label_color,
        fontsize = label_size,
        align = (:center, :center),
    )
    return nothing
end

raw"""
    _plot_points!(b::Bloch, lscene)

Plot all quantum state points on the Bloch sphere.

# Arguments
- `b::Bloch`: Contains point data and styling information
- `lscene`: LScene object for plotting

Handles both scatter points and line traces based on style specifications.
"""
function _plot_points!(b::Bloch, lscene)
    isempty(b.points) && return nothing
    for k in 1:length(b.points)
        pts = b.points[k]
        style = b.point_style[k]
        alpha = b.point_alpha[k]
        marker = b.point_marker[mod1(k, length(b.point_marker))]
        N = size(pts, 2)

        raw_x = pts[1, :]
        raw_y = pts[2, :]
        raw_z = pts[3, :]

        ds = vec(sqrt.(sum(abs2, pts; dims = 1)))
        if !all(isapprox.(ds, ds[1]; rtol = 1e-12))
            indperm = sortperm(ds)
        else
            indperm = collect(1:N)
        end
        this_color = b.point_color[k]
        if style == :m
            defaults = b.point_default_color
            L = length(defaults)
            times = ceil(Int, N / L)
            big_colors = repeat(b.point_default_color, times)[1:N]
            big_colors = big_colors[indperm]
            colors = big_colors
        else
            if this_color === nothing
                defaults = b.point_default_color
                colors = defaults[mod1(k, length(defaults))]
            else
                colors = this_color
            end
        end
        if style in (:s, :m)
            scatter!(
                lscene,
                raw_x[indperm],
                raw_y[indperm],
                raw_z[indperm];
                color = colors,
                markersize = b.point_size[mod1(k, length(b.point_size))],
                marker = marker,
                transparency = alpha < 1.0,
                alpha = alpha,
                strokewidth = 0.0,
            )

        elseif style == :l
            c = isa(colors, Vector) ? colors[1] : colors
            lines!(lscene, raw_x, raw_y, raw_z; color = c, transparency = alpha < 1.0, alpha = alpha)
        end
    end
    return nothing
end

raw"""
    _plot_lines!(b::Bloch, lscene)

Draw all connecting lines between points on the Bloch sphere.

# Arguments
- `b::Bloch`: Contains line data and formatting
- `lscene`: LScene object for drawing

Processes line style specifications and color mappings.
"""
function _plot_lines!(b::Bloch, lscene)
    isempty(b.lines) && return nothing
    color_map =
        Dict("k" => :black, "r" => :red, "g" => :green, "b" => :blue, "c" => :cyan, "m" => :magenta, "y" => :yellow)
    for (line, fmt) in b.lines
        x, y, z = line
        color_char = first(fmt)
        color = get(color_map, color_char, :black)
        linestyle = if occursin("--", fmt)
            :dash
        elseif occursin(":", fmt)
            :dot
        elseif occursin("-.", fmt)
            :dashdot
        else
            :solid
        end
        lines!(lscene, x, y, z; color = color, linestyle = linestyle)
    end
    return nothing
end

raw"""
    _plot_arcs!(b::Bloch, lscene)

Draw circular arcs connecting points on the Bloch sphere surface.

# Arguments
- `b::Bloch`: Contains arc data points
- `lscene`: LScene object for drawing

Calculates great circle arcs between specified points.
"""
function _plot_arcs!(b::Bloch, lscene)
    isempty(b.arcs) && return nothing
    for arc_pts in b.arcs
        length(arc_pts) >= 2 || continue
        v1 = normalize(arc_pts[1])
        v2 = normalize(arc_pts[end])
        n = normalize(cross(v1, v2))
        θ = acos(clamp(dot(v1, v2), -1.0, 1.0))
        if length(arc_pts) == 3
            vm = normalize(arc_pts[2])
            dot(cross(v1, vm), n) < 0 && (θ -= 2π)
        end
        t_range = range(0, θ, length = 100)
        arc_points = [Point3f((v1 * cos(t) + cross(n, v1) * sin(t))) for t in t_range]
        lines!(lscene, arc_points; color = "blue", linestyle = :solid)
    end
    return nothing
end

raw"""
    _plot_vectors!(b::Bloch, lscene)

Draw vectors from origin representing quantum states.

# Arguments
- `b::Bloch`: Contains vector data
- `lscene`: LScene object for drawing

Scales vectors appropriately and adds `3D` arrow markers.
"""
function _plot_vectors!(b::Bloch, lscene)
    isempty(b.vectors) && return nothing

    for (i, v) in enumerate(b.vectors)
        color = get(b.vector_color, i, RGBAf(0.2, 0.5, 0.8, 0.9))

        arrows3d!(
            lscene,
            [Point3f(0, 0, 0)],
            [v],
            color = color,
            shaftradius = b.vector_width,
            tiplength = b.vector_tiplength,
            tipradius = b.vector_tipradius,
            rasterize = 3,
        )
    end
    return nothing
end

@doc raw"""
    plot_bloch(::Val{:Makie}, state::QuantumObject; kwargs...)

Plot a pure quantum state on the Bloch sphere using the `Makie` backend.

# Arguments
- `state::QuantumObject{<:Union{Ket,Bra}}`: The quantum state ([`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref)) to be visualized.
- `kwargs...`: Additional keyword arguments passed to `_render_bloch_makie`.

!!! note "Internal function"
    This is the `Makie`-specific implementation called by the main `plot_bloch` function.
"""
function QuantumToolbox.plot_bloch(
    ::Val{:Makie},
    state::QuantumObject{OpType};
    kwargs...,
) where {OpType<:Union{Ket,Bra,Operator}}
    bloch_vec = _state_to_bloch(state)
    return _render_bloch_makie(bloch_vec; kwargs...)
end

end
