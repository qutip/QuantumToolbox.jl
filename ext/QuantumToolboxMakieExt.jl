module QuantumToolboxMakieExt

using QuantumToolbox
using LinearAlgebra
using Makie
using Makie:
    Axis, Axis3, Colorbar, Figure, GridLayout, heatmap!, surface!, barplot!, GridPosition, @L_str, Reverse, ylims!

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

raw"""
    _state_to_bloch(state::QuantumObject{<:Ket}) -> Vector{Float64}

Convert a pure qubit state (`Ket`) to its Bloch vector representation.

If the state is not normalized, it is automatically normalized before conversion.

# Arguments
- `state`: A `Ket` representing a pure quantum state.

# Returns
A 3-element `Vector{Float64}` representing the Bloch vector `[x, y, z]`.

# Throws
- `ArgumentError` if the state dimension is not 2.
"""
function _state_to_bloch(state::QuantumObject{<:Ket})
    if !isapprox(norm(state), 1.0, atol = 1e-6)
        @warn "State is not normalized. Normalizing before Bloch vector conversion."
        state = normalize(state)
    end
    ψ = state.data
    if length(ψ) != 2
        error("Bloch sphere visualization is only supported for qubit states (2-level systems)")
    end
    x = 2 * real(ψ[1] * conj(ψ[2]))
    y = 2 * imag(ψ[1] * conj(ψ[2]))
    z = abs2(ψ[1]) - abs2(ψ[2])
    return [x, y, z]
end

raw"""
    _dm_to_bloch(ρ::QuantumObject{<:Operator}) -> Vector{Float64}

Convert a qubit density matrix (`Operator`) to its Bloch vector representation.

This function assumes the input is Hermitian. If the density matrix is not Hermitian, a warning is issued.

# Arguments
- `ρ`: A density matrix (`Operator`) representing a mixed or pure quantum state.

# Returns
A 3-element `Vector{Float64}` representing the Bloch vector `[x, y, z]`.

# Throws
- `ArgumentError` if the matrix dimension is not 2.
"""
function _dm_to_bloch(ρ::QuantumObject{<:Operator})
    if !ishermitian(ρ)
        @warn "Density matrix is not Hermitian. Results may not be meaningful."
    end
    if size(ρ, 1) != 2
        error("Bloch sphere visualization is only supported for qubit states (2-level systems)")
    end
    x = real(ρ[1, 2] + ρ[2, 1])
    y = imag(ρ[2, 1] - ρ[1, 2])
    z = real(ρ[1, 1] - ρ[2, 2])
    return [x, y, z]
end

function _render_bloch_makie(bloch_vec::Vector{Float64}; location = nothing, kwargs...)
    b = Bloch()
    add_vectors!(b, bloch_vec)
    fig, location = _getFigAndLocation(location)
    fig, ax = render(b; location = location, kwargs...)
    return fig, ax
end

raw"""
    add_line!(
        b::QuantumToolbox.Bloch,
        start::QuantumToolbox.QuantumObject{<:Union{QuantumToolbox.Ket, QuantumToolbox.Bra, QuantumToolbox.Operator}},
        endp::QuantumToolbox.QuantumObject{<:Union{QuantumToolbox.Ket, QuantumToolbox.Bra, QuantumToolbox.Operator}};
        fmt = "k",
        kwargs...,
    )

Add a line between two quantum states or operators on the Bloch sphere visualization.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify.
- `start::QuantumObject`: The starting quantum state or operator. Can be a `Ket`, `Bra`, or `Operator`.
- `endp::QuantumObject`: The ending quantum state or operator. Can be a `Ket`, `Bra`, or `Operator`.
- `fmt::String="k"`: (optional) A format string specifying the line style and color (default is black `"k"`).
- `kwargs...`: Additional keyword arguments forwarded to the underlying line drawing function.

# Description

This function converts the given quantum objects (states or operators) into their Bloch vector representations and adds a line between these two points on the Bloch sphere visualization. 

# Example

```julia
b = Bloch()
ψ₁ = normalize(basis(2, 0) + basis(2, 1))
ψ₂ = normalize(basis(2, 0) - im * basis(2, 1))
add_line!(b, ψ₁, ψ₂; fmt = "r--")
```
"""
function QuantumToolbox.add_line!(
    b::QuantumToolbox.Bloch,
    start::QuantumObject{<:Union{Ket,Bra,Operator}},
    endp::QuantumObject{<:Union{Ket,Bra,Operator}};
    fmt = "k",
    kwargs...,
)
    p1 = if isket(start) || isbra(start)
        _state_to_bloch(start)
    else
        _dm_to_bloch(start)
    end
    p2 = if isket(endp) || isbra(endp)
        _state_to_bloch(endp)
    else
        _dm_to_bloch(endp)
    end
    return add_line!(b, p1, p2; fmt = fmt, kwargs...)
end

raw"""
    QuantumToolbox.add_states!(b::Bloch, states::QuantumObject...)

Add one or more quantum states to the Bloch sphere visualization by converting them into Bloch vectors.

# Arguments
- `b::Bloch`: The Bloch sphere object to modify
- `states::QuantumObject...`: One or more quantum states (Ket, Bra, or Operator)

"""
function QuantumToolbox.add_states!(b::Bloch, states::Vector{<:QuantumObject})
vecs = map(states) do state
    if isket(state) || isbra(state)
        return _state_to_bloch(state)
    else
        return _dm_to_bloch(state)
    end
end
    return add_vectors!(b, vecs)
end

raw"""
    render(b::QuantumToolbox.Bloch; location=nothing)

Render the Bloch sphere visualization from the given `Bloch` object `b`.

# Arguments

- `b::QuantumToolbox.Bloch`
  The Bloch sphere object containing states, vectors, and settings to visualize.

- `location` (optional)  
  Specifies where to display or save the rendered figure.
  - If `nothing` (default), the figure is displayed interactively.
  - If a file path (String), the figure is saved to the specified location.
  - Other values depend on backend support.

# Returns

- A tuple `(fig, axis)` where `fig` is the figure object and `axis` is the axis object used for plotting.
  These can be further manipulated or saved by the user.
"""
function QuantumToolbox.render(b::QuantumToolbox.Bloch; location = nothing)
    fig, location = _getFigAndLocation(location)
    if isnothing(b.fig) || b.fig !== fig
        b.fig = fig
        b.ax = Makie.Axis3(
            location;
            aspect = :data,
            limits = (-1.1, 1.1, -1.1, 1.1, -1.1, 1.1),
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
            zticklabelsvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
            xlabel = "",
            ylabel = "",
            zlabel = "",
            backgroundcolor = Makie.RGBAf(1, 1, 1, 0.0),
            xypanelvisible = false,
            xzpanelvisible = false,
            yzpanelvisible = false,
            xspinesvisible = false,
            yspinesvisible = false,
            zspinesvisible = false,
            protrusions = (0, 0, 0, 0),
            viewmode = :fit,
        )
        b.ax.azimuth[] = deg2rad(b.view_angles[1])
        b.ax.elevation[] = deg2rad(b.view_angles[2])
    else
        if !(b.ax in contents(location))
            location[] = b.ax
        end
    end
    Makie.empty!(b.ax)
    sphere_color = Makie.RGBAf(1.0, 0.86, 0.86, 0.2)
    Makie.mesh!(b.ax, Makie.Sphere(Point3f(0, 0, 0), 1.0f0); color = sphere_color, transparency = false)
    wire_color = Makie.RGBAf(0.5, 0.5, 0.5, 0.4)
    φ = range(0, 2π, length = 100)
    Makie.lines!(b.ax, [Point3f(sin(φi), -cos(φi), 0) for φi in φ]; color = wire_color, linewidth = 1.0)
    Makie.lines!(b.ax, [Point3f(0, -cos(φi), sin(φi)) for φi in φ]; color = wire_color, linewidth = 1.0)
    Makie.lines!(b.ax, [Point3f(sin(φi), 0, cos(φi)) for φi in φ]; color = wire_color, linewidth = 1.0)
    axis_color = Makie.RGBAf(0.3, 0.3, 0.3, 0.8)
    axis_width = 0.8
    Makie.lines!(b.ax, [Point3f(0, -1.01, 0), Point3f(0, 1.01, 0)]; color = axis_color, linewidth = axis_width)
    Makie.lines!(b.ax, [Point3f(-1.01, 0, 0), Point3f(1.01, 0, 0)]; color = axis_color, linewidth = axis_width)
    Makie.lines!(b.ax, [Point3f(0, 0, -1.01), Point3f(0, 0, 1.01)]; color = axis_color, linewidth = axis_width)
    for (k, points) in enumerate(b.points)
        style = b.point_style[k]
        color = b.point_color[k]
        alpha = b.point_alpha[k]
        marker = b.point_marker[(k-1)%length(b.point_marker)+1]
        x = points[2, :]
        y = -points[1, :]
        z = points[3, :]

        if style == :s || style == :m
            Makie.scatter!(
                b.ax,
                x,
                y,
                z;
                color = color,
                markersize = 6,
                marker = marker,
                strokewidth = 0.05,
                transparency = alpha < 1,
                alpha = alpha,
            )
        elseif style == :l
            Makie.lines!(b.ax, x, y, z; color = color, linewidth = 2.0, transparency = alpha < 1, alpha = alpha)
        end
    end
    for (line, fmt, kwargs) in b.lines
        x, y, z = line
        color_map =
            Dict("k" => :black, "r" => :red, "g" => :green, "b" => :blue, "c" => :cyan, "m" => :magenta, "y" => :yellow)
        c = get(color_map, first(fmt), :black)
        ls = nothing
        if occursin("--", fmt)
            ls = :dash
        elseif occursin(":", fmt)
            ls = :dot
        elseif occursin("-.", fmt)
            ls = :dashdot
        else
            ls = :solid
        end
        color = get(kwargs, :color, c)
        linewidth = get(kwargs, :linewidth, 1.0)
        linestyle = get(kwargs, :linestyle, ls)
        Makie.lines!(b.ax, x, y, z; color = color, linewidth = linewidth, linestyle = linestyle)
    end
    for arc_pts in b.arcs
        if length(arc_pts) >= 2
            v1 = normalize(arc_pts[1])
            v2 = normalize(arc_pts[end])
            n = normalize(cross(v1, v2))
            θ = acos(clamp(dot(v1, v2), -1.0, 1.0))
            if length(arc_pts) == 3
                vm = normalize(arc_pts[2])
                if dot(cross(v1, vm), n) < 0
                    θ = 2π - θ
                end
            end
            t_range = range(0, θ, length = 100)
            arc_points = [Makie.Point3f((v1*cos(t) + cross(n, v1)*sin(t))...) for t in t_range]
            Makie.lines!(b.ax, arc_points; color = Makie.RGBAf(0.8, 0.4, 0.1, 0.9), linewidth = 2.0, linestyle = :solid)
        end
    end
    r = 1.0
    if !isempty(b.vectors)
        for (i, v) in enumerate(b.vectors)
            color = Base.get(b.vector_color, i, Makie.RGBAf(0.2, 0.5, 0.8, 0.9))
            start = Point3f(0, 0, 0)
            vec = Vec3f(v[2], -v[1], v[3])
            length = norm(vec)
            max_length = r * 0.79
            if length > max_length
                vec = (vec / length) * max_length
            end
            endp = Point3f(vec)
            arrowsize = Vec3f(0.07, 0.08, 0.08)
            Makie.arrows!(
                b.ax,
                [start],
                [vec];
                color = color,
                linewidth = 0.028,
                arrowsize = arrowsize,
                arrowcolor = color,
            )
        end
    end
    label_color = Makie.RGBf(0.2, 0.2, 0.2)
    label_size = 16
    label_font = "TeX Gyre Heros"
    Makie.text!(
        b.ax,
        "y";
        position = Point3f(1.04, 0, 0),
        color = label_color,
        fontsize = label_size,
        font = label_font,
    )
    Makie.text!(
        b.ax,
        "x";
        position = Point3f(0, -1.10, 0),
        color = label_color,
        fontsize = label_size,
        font = label_font,
    )
    Makie.text!(
        b.ax,
        "|1⟩";
        position = Point3f(0, 0, -1.10),
        color = label_color,
        fontsize = label_size,
        font = label_font,
    )
    Makie.text!(
        b.ax,
        "|0⟩";
        position = Point3f(0, 0, 1.08),
        color = label_color,
        fontsize = label_size,
        font = label_font,
    )
    return fig, b.ax
end

function QuantumToolbox.plot_bloch(::Val{:Makie}, state::QuantumObject{<:Union{Ket,Bra}}; kwargs...)
    state = isbra(state) ? dag(state) : state
    bloch_vec = _state_to_bloch(state)
    return _render_bloch_makie(bloch_vec; kwargs...)
end

function QuantumToolbox.plot_bloch(::Val{:Makie}, ρ::QuantumObject{<:Operator}; kwargs...)
    bloch_vec = _dm_to_bloch(ρ)
    return _render_bloch_makie(bloch_vec; kwargs...)
end

end
