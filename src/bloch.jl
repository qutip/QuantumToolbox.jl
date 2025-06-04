export Bloch, set_label_convention!, clear!, add_points!, add_states!, add_vectors!, add_annotation!, add_line!, add_arc!, render!, show!, save

_state_to_cartesian_coordinates(state::Qobj) = [
    real(expect(sigmax(), state)),
    real(expect(sigmay(), state)),
    real(expect(sigmaz(), state))
]

@doc raw"""
    Bloch(fig=nothing, axes=nothing, view=nothing, figsize=nothing, background=false)

Class for plotting data on the Bloch sphere.  Valid data can be either
points, vectors, or Qobj objects.

Attributes
----------
fig : Makie.Figure
    Makie Figure instance for plotting Bloch sphere.
axes : Makie.LScene
    Makie LScene for Bloch sphere animation.
font_color : default "black"
    Color of font used for Bloch sphere labels.
font_size : Int, default 20
    Size of font used for Bloch sphere labels.
frame_alpha : Float64, default 0.2
    Sets transparency of Bloch sphere frame.
frame_color : default "gray"
    Color of sphere wireframe.
frame_width : Int, default 1
    Width of wireframe.
point_color : Vector, default ["blue", "red", "green", "#CC6600"]
    List of colors for Bloch sphere point markers to cycle through, i.e.
    by default, points 1 and 5 will both be blue ("blue").
point_marker : Vector, default [:circle, :rect, :diamond, :utriangle]
    List of point marker shapes to cycle through.
point_size : Vector{Int}, default [12, 14, 12, 14]
    List of point marker sizes.
sphere_alpha : Float64, default 0.15
    Transparency of Bloch sphere itself.
sphere_color : default "#FFDDDD"
    Color of Bloch sphere.
figsize : Tuple, default [5, 5]
    Figure size of Bloch sphere plot.  Best to have both numbers the same;
vector_color : Vector, ["green", "#CC6600", "blue", "red"]
    List of vector colors to cycle through.
vector_width : Float64, default 0.03
    Width of displayed vectors.
vector_style : Char, default '▲'
    Vector arrowhead style
vector_mutation : Float64, default 0.1
    Width of vectors arrowhead.
view : Vector, default [-60, 30]
    Azimuthal and Elevation viewing angles.
xlabel : Vector, default ["x", ""]
    List of strings corresponding to +x and -x axes labels, respectively.
xlpos : Vector, default [1.2, -1.2]
    Positions of +x and -x labels respectively.
ylabel : Vector, default ["y", ""]
    List of strings corresponding to +y and -y axes labels, respectively.
ylpos : Vector, default [1.2, -1.2]
    Positions of +y and -y labels respectively.
zlabel : Vector, default ["|0⟩", "|1⟩"]
    List of strings corresponding to +z and -z axes labels, respectively.
zlpos : Vector, default [1.2, -1.2]
    Positions of +z and -z labels respectively.
"""
mutable struct Bloch
    # Figure and axes
    fig
    axes
    # Background axes, default = false
    background::Bool
    # The size of the figure in inches, default = [5,5].
    figsize::Tuple
    # Azimuthal and Elvation viewing angles, default = [-60,30].
    view::Vector
    # Color of Bloch sphere, default = #FFDDDD
    sphere_color
    # Transparency of Bloch sphere, default = 0.2
    sphere_alpha::Float64
    # Color of wireframe, default = 'gray'
    frame_color
    # Width of wireframe, default = 1
    frame_width::Float64
    # Transparency of wireframe, default = 0.2
    frame_alpha::Float64
    # Labels for x-axis (in LaTex), default = ['$x$', '']
    xlabel::Vector
    # Position of x-axis labels, default = [1.2, -1.2]
    xlpos::Vector{Float64}
    # Labels for y-axis (in LaTex), default = ['$y$', '']
    ylabel::Vector
    # Position of y-axis labels, default = [1.2, -1.2]
    ylpos::Vector{Float64}
    # Labels for z-axis (in LaTex),
    # default = [r'$\left\|0\right>$', r'$\left|1\right>$']
    zlabel::Vector
    # Position of z-axis labels, default = [1.2, -1.2]
    zlpos::Vector{Float64}
    # ---font options---
    # Color of fonts, default = 'black'
    font_color
    # Size of fonts, default = 20
    font_size::Int

    # ---vector options---
    # List of colors for Bloch vectors, default = ['b','g','r','y']
    vector_default_color::Vector
    # List that stores the display colors for each vector
    vector_color::Vector
    # Width of Bloch vectors, default = 5
    vector_width::Float64
    # Style of Bloch vectors, default = '-\|>' (or 'simple')
    vector_style::Char
    # Sets the width of the vectors arrowhead
    vector_mutation::Float64

    # ---point options---
    # List of colors for Bloch point markers, default = ['b','g','r','y']
    point_default_color::Vector
    # Old variable used in V4 to customise the color of the points
    point_color
    # List that stores the display colors for each set of points
    _inner_point_color::Vector
    # Size of point markers, default = 25
    point_size::Vector
    # Shape of point markers, default = ['o','^','d','s']
    point_marker::Vector

    # ---data lists---
    # Data for point markers
    points::Vector
    # Data for Bloch vectors
    vectors::Vector
    # Transparency of vectors, alpha value from 0 to 1
    vector_alpha::Vector{Float64}
    # Data for annotations
    annotations::Vector
    # Number of times sphere has been saved
    savenum::Int
    # Style of points, 'm' for multiple colors, 's' for single color
    point_style::Vector
    # Transparency of points, alpha value from 0 to 1
    point_alpha::Vector{Float64}
    # Data for line segment
    _lines::Vector
    # Data for arcs and arc style
    _arcs::Vector

    function Bloch(; fig=nothing, axes=nothing, view=nothing, figsize=nothing, background=false)
        new(
            fig,                                    # fig    
            axes,                                   # axes
            background,                             # background
            isnothing(figsize) ? (5, 5) : figsize,  # figsize
            isnothing(view) ? [-60, 30] : view,     # view
            "#FFDDDD",                              # sphere_color
            0.15,                                   # sphere_alpha
            "gray",                                 # frame_color
            1,                                      # frame_width
            0.2,                                    # frame_alpha
            ["x", ""],                              # xlabel
            [1.2, -1.2],                            # xlpos
            ["y", ""],                              # ylabel
            [1.2, -1.2],                            # ylpos
            ["|0⟩", "|1⟩"],                          # zlabel
            [1.2, -1.2],                            # zlpos
            "black",                                # font_color
            20,                                     # font_size
            ["green", "#CC6600", "blue", "red"],    # vector_default_color
            [],                                     # vector_color
            0.03,                                   # vector_width
            '▲',                                    # vector_style
            0.1,                                    # vector_mutation
            ["blue", "red", "green", "#CC6600"],    # point_default_color
            nothing,                                # point_color
            [],                                     # _inner_point_color
            [12, 14, 12, 14],                       # point_size
            [:circle, :rect, :diamond, :utriangle], # point_marker
            [],                                     # points
            [],                                     # vectors
            Float64[],                              # vector_alpha
            [],                                     # annotations
            0,                                      # savenum
            [],                                     # point_style
            Float64[],                              # point_alpha
            [],                                     # _lines
            []                                      # _arcs
        )
    end
end

@doc raw"""
    set_label_convention!(b::Bloch, convention::String)

Set x, y and z labels according to one of conventions.

Parameters
----------
convention : String
    One of the following:

    - "original"
    - "xyz"
    - "sx sy sz"
    - "01"
    - "polarization jones"
    - "polarization jones letters"
        see also: https://en.wikipedia.org/wiki/Jones_calculus
    - "polarization stokes"
        see also: https://en.wikipedia.org/wiki/Stokes_parameters
"""
function set_label_convention!(b::Bloch, convention::String)
    if convention == "original"
        b.xlabel = ["x", ""]
        b.ylabel = ["y", ""]
        b.zlabel = ["|0⟩", "|1⟩"]
    elseif convention == "xyz"
        b.xlabel = ["x", ""]
        b.ylabel = ["y", ""]
        b.zlabel = ["z", ""]
    elseif convention == "sx sy sz"
        b.xlabel = [L"s_x", ""]
        b.ylabel = [L"s_y", ""]
        b.zlabel = [L"s_z", ""]
    elseif convention == "01"
        b.xlabel = ["", ""]
        b.ylabel = ["", ""]
        b.zlabel = ["|0⟩", "|1⟩"]
    elseif convention == "polarization jones" 
        b.xlabel = ["|↗︎↙︎⟩", "|↖︎↘︎⟩"]     # TODO: fix spacing
        b.ylabel = ["|↺⟩", "|↻⟩"]
        b.zlabel = ["|↔︎⟩", "|↕︎⟩"]
    elseif convention == "polarization jones letters"
        b.xlabel = ["|D⟩", "|A⟩"]
        b.ylabel = ["|L⟩", "|R⟩"]
        b.zlabel = ["|H⟩", "|V⟩"]
    elseif convention == "polarization stokes"
        b.ylabel = ["↗︎↙︎", "↖︎↘︎"]     # TODO: fix spacing
        b.zlabel = ["↺", "↻"]
        b.xlabel = ["↔︎", "↕︎"]
    else
        error("No such convention.")
    end
end

function Base.show(io::IO, b::Bloch)
    println(io, "Bloch data:")
    println(io, "-----------")
    println(io, "Number of points:  ", length(b.points))
    println(io, "Number of vectors: ", length(b.vectors))
    println(io)
    println(io, "Bloch sphere properties:")
    println(io, "------------------------")
    println(io, "font_color:           ", b.font_color)
    println(io, "font_size:            ", b.font_size)
    println(io, "frame_alpha:          ", b.frame_alpha)
    println(io, "frame_color:          ", b.frame_color)
    println(io, "frame_width:          ", b.frame_width)
    println(io, "point_default_color:  ", b.point_default_color)
    println(io, "point_marker:         ", b.point_marker)
    println(io, "point_size:           ", b.point_size)
    println(io, "sphere_alpha:         ", b.sphere_alpha)
    println(io, "sphere_color:         ", b.sphere_color)
    println(io, "figsize:              ", b.figsize)
    println(io, "vector_default_color: ", b.vector_default_color)
    println(io, "vector_width:         ", b.vector_width)
    println(io, "vector_style:         ", b.vector_style)
    println(io, "vector_mutation:      ", b.vector_mutation)
    println(io, "view:                 ", b.view)
    println(io, "xlabel:               ", b.xlabel)
    println(io, "xlpos:                ", b.xlpos)
    println(io, "ylabel:               ", b.ylabel)
    println(io, "ylpos:                ", b.ylpos)
    println(io, "zlabel:               ", b.zlabel)
    println(io, "zlpos:                ", b.zlpos)
end
function Base.show(io::IO, ::MIME"image/png", b::Bloch)
    render!(b)
    buf = IOBuffer()
    save(buf, MIME("image/png"), b.fig)
    write(io, take!(buf))
end
function Base.show(io::IO, ::MIME"image/svg+xml", b::Bloch)
    render!(b)
    buf = IOBuffer()
    save(buf, MIME("image/svg+xml"), b.fig)
    write(io, String(take!(buf)))
end

@doc raw"""
    clear!(b::Bloch)    

Resets Bloch sphere data sets to empty.
"""
function clear!(b::Bloch)
    empty!(b.points)
    empty!(b.vectors)
    empty!(b.point_style)
    empty!(b.point_alpha)
    empty!(b.vector_alpha)
    empty!(b.annotations)
    empty!(b.vector_color)
    b.point_color = nothing
    empty!(b._lines)
    empty!(b._arcs)
end

@doc raw"""
    add_points!(b::Bloch, points; meth="s", colors=nothing, alpha::Float64=1.0)

Add a list of data points to bloch sphere.

Parameters
----------
points : Array
    Collection of data points.

meth : {'s', 'm', 'l'}
    Type of points to plot, use 'm' for multicolored, 'l' for points
    connected with a line.

colors : 
    Optional array with colors for the points.
    A single color for meth 's', and list of colors for meth 'm'

alpha : Float64, default=1.0
    Transparency value for the vectors. Values between 0 and 1.
"""
function add_points!(b::Bloch, points::Array; meth="s", colors=nothing, alpha::Float64=1.0)
    if ndims(points) == 1
        if eltype(points) <: Number
            points = reshape(points, 3,1)
        elseif eltype(points) <: AbstractVector
            points = reduce(hcat, points)
        end
    end

    if ndims(points) != 2 || size(points, 1) != 3
        error("The included points are not valid. Points must be equivalent to a 2D array where the first index represents the x,y,z values and the second index iterates over the points.")
    end

    if meth ∉ ["s", "m", "l"]
        error("The value for meth = $meth is not valid. Please use 's', 'l' or 'm'.")
    end

    push!(b.point_style, meth)
    push!(b.points, points)
    push!(b.point_alpha, alpha)
    push!(b._inner_point_color, colors)
end

@doc raw"""
    add_states!(b::Bloch, state::Qobj; kind="vector", colors=nothing, alpha=1.0)

Add a state vector Qobj to Bloch sphere.

Parameters
----------
state : Vector{Qobj} or Qobj
    Input state vector or list.

kind : {'vector', 'point'}
    Type of object to plot.

colors :
    Optional array with colors for the states.
    The colors can be a string or a RGB or RGBA tuple.

alpha : Float64, default=1.0
    Transparency value for the vectors. Values between 0 and 1.
"""
function add_states!(b::Bloch, state::Qobj; kind="vector", colors=nothing, alpha=1.0)
    vec = _state_to_cartesian_coordinates(state)
    if kind == "vector"
        add_vectors!(b, vec, colors=colors, alpha=alpha)
    elseif kind == "point"
        add_points!(b, vec, colors=colors, alpha=alpha)
    else
        error("The included kind is not valid. It should be vector or point, not $kind.")
    end
end
function add_states!(b::Bloch, states::Vector{<:Qobj}; kind="vector", colors=nothing, alpha=1.0)
    if !isnothing(colors)
        if ndims(colors) == 0
            colors = fill(colors, length(states))
        elseif ndims(colors) == 1 && eltype(colors) <: Number
            colors = fill(colors, length(states))
        end

        if length(colors) != length(states) || ndims(colors) > 2 || ndims(colors) == 2 && !(eltype(colors) <: Number)
            error("The included colors are not valid. colors must have the same size as state.")
        end
    else
        colors = fill(nothing, length(states))
    end

    for (k, st) in enumerate(states)
        add_states!(b, st; kind=kind, colors=colors[k], alpha=alpha)
    end
end

@doc raw"""
    add_vectors!(b::Bloch, vectors; colors=nothing, alpha=1.0)

Add a list of vectors to Bloch sphere.

Parameters
----------
vectors :
    Array with vectors of unit length or smaller.

colors :
    Optional array with colors for the vectors.
    The colors can be a string or a RGB or RGBA tuple.

alpha : Float64, default=1.0
    Transparency value for the vectors. Values between 0 and 1.
"""
function add_vectors!(b::Bloch, vectors; colors=nothing, alpha::Float64=1.0)
    if ndims(vectors) == 1
        if eltype(vectors) <: Number
            vectors = reshape(vectors, 1,3)
        elseif eltype(vectors) <: AbstractVector
            vectors = reduce(hcat, vectors)'
        end
    end

    if ndims(vectors) != 2 || size(vectors, 2) != 3
        error("The included vectors are not valid. Vectors must be equivalent to a 2D array where the first index represents the iteration over the vectors and the second index represents the position in 3D of vector head.")
    end

    colors = nothing
    if isnothing(colors)
        colors = fill(nothing, size(vectors, 1))
    else
        if ndims(colors) == 0
            colors = fill(colors, size(vectors, 1))
        end

        if size(colors, 1) != size(vectors, 1) || ndims(colors) == 2 && !(eltype(colors) <: Number)
            error("The included colors are not valid. colors must have the same size as vectors.")
        end
    end

    for (k, vec) in enumerate(eachrow(vectors))
        push!(b.vectors, vec)
        push!(b.vector_alpha, alpha)
        push!(b.vector_color, colors[k])
    end
end

@doc raw"""
    add_annotation!(b::Bloch, state_or_vector, text; kwargs...)

Add a text or LaTeX annotation to Bloch sphere, parametrized by a qubit
state or a vector.

Parameters
----------
state_or_vec : Qobj or Vector
    Position for the annotaion.
    Qobj of a qubit or a vector of 3 elements.

text : String
    Annotation text.
    You can use LaTeX, but remember to use L string
    e.g. L"\langle x \rangle"

kwargs :
    Options as for Makie.text!, including:
    fontsize, color, align.
"""
function add_annotation!(b::Bloch, state::Qobj, text::String; kwargs...)
    vec = _state_to_cartesian_coordinates(state)
    add_annotation!(b, vec, text; kwargs...)
end
function add_annotation!(b::Bloch, vec, text::String; kwargs...)
    if length(vec) != 3
        error("Position needs to be specified by a qubit state or a 3D vector.")
    end

    push!(b.annotations, Dict(
        :position => Point3f(vec),
        :text => text,
        :opts => kwargs
    ))
end

@doc raw"""
    add_arc!(b::Bloch, start, endpt; fmt=:solid, steps=nothing, kwargs...)

Adds an arc between two points on a sphere. The arc is set to be
blue solid curve by default.

The start and end points must be on the same sphere (i.e. have the
same radius) but need not be on the unit sphere.

Parameters
----------
start : Qobj or Vector
    Array with cartesian coordinates of the first point, or a state
    vector or density matrix that can be mapped to a point on or
    within the Bloch sphere.
end : Qobj or Vector
    Array with cartesian coordinates of the second point, or a state
    vector or density matrix that can be mapped to a point on or
    within the Bloch sphere.
fmt : Symbol, default: :solid
    A Makie format symbol for rendering the arc.
steps : Int, default: None
    The number of segments to use when rendering the arc. The default
    uses 100 steps times the distance between the start and end points,
    with a minimum of 2 steps.
kwargs... : Dict
    Additional parameters to pass to the Makie.lines! function
    when rendering this arc.
"""
function add_arc!(b::Bloch, start, endpt; fmt=:solid, steps=nothing, kwargs...)
    pt1 = start isa Qobj ? _state_to_cartesian_coordinates(start) : start
    pt2 = endpt isa Qobj ? _state_to_cartesian_coordinates(endpt) : endpt

    len1 = norm(pt1)
    len2 = norm(pt2)

    if len1 < 1e-12 || len2 < 1e-12
        error("Polar and azimuthal angles undefined at origin.")
    elseif abs(len1 - len2) > 1e-12
        error("Points not on the same sphere.")
    elseif norm(pt1 - pt2) < 1e-12
        error("Start and end represent the same point. No arc can be formed.")
    elseif norm(pt1 + pt2) < 1e-12
        error("Start and end are diagonally opposite, no unique arc is possible.")
    end
    
    if isnothing(steps)
        steps = round(Int, norm(pt1 - pt2) * 100)
        steps = max(2, steps)
    end
    t = range(0, 1; length=steps)
    # All the points in this line are contained in the plane defined
    # by pt1, pt2 and the origin.
    line = [pt1 .* τ + pt2 .* (1-τ) for τ in t]

    # Normalize all the points in the line so that are distance len1 from
    # the origin.
    arc = Point3f.([v .* len1 / norm(v) for v in line])

    # Save the arc
    push!(b._arcs, (arc, fmt, Dict(kwargs)))
end

@doc raw"""
    add_line!(b::Bloch, start, endpt, fmt=:solid, kwargs...)

Adds a line segment connecting two points on the bloch sphere.

The line segment is set to be a black solid line by default.

Parameters
----------
start : Qobj or Vector
    Array with cartesian coordinates of the first point, or a state
    vector or density matrix that can be mapped to a point on or
    within the Bloch sphere.
end : Qobj or Vector
    Array with cartesian coordinates of the second point, or a state
    vector or density matrix that can be mapped to a point on or
    within the Bloch sphere.
fmt : Symbol, default: :solid
    A Makie format symbol for rendering the line.
**kwargs : Dict
    Additional parameters to pass to the matplotlib .plot function
    when rendering this line.
"""
function add_line!(b::Bloch, start, endpt; fmt=:solid, kwargs...)
    pt1 = start isa Qobj ? _state_to_cartesian_coordinates(start) : start
    pt2 = endpt isa Qobj ? _state_to_cartesian_coordinates(endpt) : endpt

    push!(b._lines, ([Point3f(pt1), Point3f(pt2)], fmt, kwargs))
end

@doc raw"""
    render!(b::Bloch)

Render the Bloch sphere and its data sets in on given figure and axes.
"""
function render!(b::Bloch)
    if isnothing(b.fig)
        b.fig = Figure(size=b.figsize.*96)
    end

    if isnothing(b.axes)
        b.axes = LScene(b.fig[1,1], show_axis=false)
        update_cam!(b.axes.scene, (0,1,0), (0,0,0))
        rotate_cam!(b.axes.scene, (-deg2rad(b.view[2]), deg2rad(b.view[1]), 0))
    end

    plot_sphere!(b)
    plot_points!(b)
    plot_vectors!(b)
    plot_lines!(b)
    plot_arcs!(b)
    if !b.background
        plot_axes!(b)
    end
    plot_axes_labels!(b)
    plot_annotations!(b)
end

function plot_sphere!(b::Bloch)
    u = range(0, 2π, length=50)
    v = range(0, π, length=25)
    x = cos.(u) * transpose(sin.(v))
    y = sin.(u) * transpose(sin.(v))
    z = ones(size(u)) * transpose(cos.(v))

    # sphere
    mesh!(b.axes, Sphere(Point3f(0), 1f0), color=b.sphere_color, transparency=true, alpha=b.sphere_alpha)

    # wireframe
    for i in 1:5:50
        lines!(b.axes, x[i, :], y[i, :], z[i, :], color=b.frame_color, transparency=true, alpha=b.frame_alpha, linewidth=b.frame_width)
    end
    for j in 1:5:25
        lines!(b.axes, x[:, j], y[:, j], z[:, j], color=b.frame_color, transparency=true, alpha=b.frame_alpha, linewidth=b.frame_width)
    end

    # equator
    lines!(b.axes, cos.(u), sin.(u), zeros(size(u)), color=b.frame_color, linewidth=b.frame_width, transparency=true)
    lines!(b.axes, zeros(size(u)), cos.(u), sin.(u), color=b.frame_color, linewidth=b.frame_width, transparency=true)
end

function plot_axes!(b::Bloch)
    # axes
    lines!(b.axes, [Point3f(1,0,0), Point3f(-1,0,0)], color=b.frame_color, linewidth=b.frame_width, transparency=true)
    lines!(b.axes, [Point3f(0,1,0), Point3f(0,-1,0)], color=b.frame_color, linewidth=b.frame_width, transparency=true)
    lines!(b.axes, [Point3f(0,0,1), Point3f(0,0,-1)], color=b.frame_color, linewidth=b.frame_width, transparency=true)
end

function plot_axes_labels!(b::Bloch)
    # axes labels
    text!(b.axes, b.xlabel[1], position = Vec3f(b.xlpos[1],0,0), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)
    text!(b.axes, b.xlabel[2], position = Vec3f(b.xlpos[2],0,0), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)

    text!(b.axes, b.ylabel[1], position = Vec3f(0,b.ylpos[1],0), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)
    text!(b.axes, b.ylabel[2], position = Vec3f(0,b.ylpos[2],0), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)

    text!(b.axes, b.zlabel[1], position = Vec3f(0,0,b.zlpos[1]), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)
    text!(b.axes, b.zlabel[2], position = Vec3f(0,0,b.zlpos[2]), align=(:center,:center), color=b.font_color, fontsize=b.font_size, transparency=true)
end

function plot_vectors!(b::Bloch)  # TODO: alpha not working
    for (k, vec) in enumerate(b.vectors)
        alpha = b.vector_alpha[k]
        color = b.vector_color[k]
        if isnothing(color)
            idx = mod1(k, length(b.vector_default_color))
            color = b.vector_default_color[idx]
        end

        arrows!(b.axes, [Point3f(0)], [Point3f(vec)], linewidth=b.vector_width, color=color, transparency=true, alpha=alpha, arrowsize=b.vector_mutation, fxaa=true)
    end
end

function plot_points!(b::Bloch)
    for (k, points) in enumerate(b.points)
        num_points = size(points, 2)

        dist = [norm(points[:,j]) for j in 1:num_points]
        if !all(isapprox.(dist, dist[1]; rtol=1e-12))
            indperm = sortperm(dist)
        else
            indperm = 1:num_points
        end

        s = b.point_size[mod1(k, length(b.point_size))]
        marker = b.point_marker[mod1(k, length(b.point_marker))]
        style = b.point_style[k]

        color = nothing
        if !isnothing(b._inner_point_color[k])
            color = self._inner_point_color[k]
        elseif !isnothing(b.point_color)
            color = b.point_color
        elseif style ∈ ["s", "l"]
            color = b.point_default_color[mod1(k, length(b.point_default_color))]
        elseif style == "m"
            nrepeats = ceil(Int, num_points / length(b.point_default_color))
            color = repeat(b.point_default_color, nrepeats)[1:num_points]
            color = color[indperm]
        end

        if style ∈ ["s", "m"]
            scatter!(b.axes, points[1, indperm], points[2, indperm], points[3, indperm], color=color, markersize=s, marker=marker, transparency=true, alpha=b.point_alpha[k])
        elseif style == "l"
            lines!(b.axes, points[1, :], points[2, :], points[3, :], color=color, transparency=true, alpha=b.point_alpha[k])
        end
    end
end

function plot_annotations!(b::Bloch)
    for annotation in b.annotations
        opts = Dict(:fontsize => b.font_size, :color => b.font_color, :transparency => true)
        merge!(opts, annotation[:opts])
        text!(b.axes, annotation[:text], position=annotation[:position], align=(:center, :center); opts...)
    end
end

function plot_lines!(b::Bloch)
    for (line, fmt, kw) in b._lines
        lines!(b.axes, line; linestyle=fmt, transparency=true, kw...)
    end
end

function plot_arcs!(b::Bloch)
    for (arc, fmt, kw) in b._arcs
        lines!(b.axes, arc; linestyle=fmt, transparency=true, kw...)
    end
end

@doc raw"""
    show!(b::Bloch)

Display Bloch sphere and corresponding data sets.
"""
function show!(b::Bloch)
    render!(b)
    display(b.fig)
end

@doc raw"""
    save(b::Bloch; name=nothing, format="png", dirc=nothing, dpin=nothing)

Saves Bloch sphere to file of type ``format`` in directory ``dirc``.

Parameters
----------

name : String
    Name of saved image. Must include path and format as well.
    i.e. '/Users/Me/Desktop/bloch.png'
    This overrides the 'format' and 'dirc' arguments.
format : String
    Format of output image.
dirc : String
    Directory for output images. Defaults to current working directory.
"""
function save(b::Bloch, name::String=nothing, format::String="png", dirc::String=nothing)
    render!(b)
    # Conditional variable for first argument to savefig
    # that is set in subsequent if-elses
    if !(isnothing(dirc))
        dirpath = joinpath(pwd(), dirc)
        if !isdir(dirpath)
            mkpath(dirpath)
        end
    else
        dirpath = pwd()
    end
    if isnothing(name)
        complete_path = joinpath(dirpath, "bloch_$(b.savenum).$format")
    else
        complete_path = name
    end

    save(complete_path, b.fig)
    b.savenum += 1
end