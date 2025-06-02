export plot_wigner,
    plot_fock_distribution,
    plot_bloch,
    Bloch,
    add_points!,
    add_vectors!,
    add_line!,
    add_arc!,
    clear!,
    render,
    add_states!

@doc raw"""
    plot_wigner(
        state::QuantumObject{OpType}; 
        library::Union{Val,Symbol}=Val(:Makie), 
        kwargs...
    ) where {OpType<:Union{Bra,Ket,Operator}

Plot the [Wigner quasipropability distribution](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) of `state` using the [`wigner`](@ref) function. 
    
The `library` keyword argument specifies the plotting library to use, defaulting to [`Makie.jl`](https://github.com/MakieOrg/Makie.jl). 

# Arguments
- `state::QuantumObject`: The quantum state for which to plot the Wigner distribution.
- `library::Union{Val,Symbol}`: The plotting library to use. Default is `Val(:Makie)`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. See the documentation for the specific plotting library for more information.

!!! note "Import library first"
    The plotting libraries must first be imported before using them with this function.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:Makie)` instead of `:Makie` as the plotting library. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
plot_wigner(
    state::QuantumObject{OpType};
    library::Union{Val,Symbol} = Val(:Makie),
    kwargs...,
) where {OpType<:Union{Bra,Ket,Operator}} = plot_wigner(makeVal(library), state; kwargs...)

plot_wigner(::Val{T}, state::QuantumObject{OpType}; kwargs...) where {T,OpType<:Union{Bra,Ket,Operator}} =
    throw(ArgumentError("The specified plotting library $T is not available. Try running `using $T` first."))

@doc raw"""
    plot_fock_distribution(
        ρ::QuantumObject{SType};
        library::Union{Val, Symbol} = Val(:Makie),
        kwargs...
    ) where {SType<:Union{Ket,Operator}}

Plot the [Fock state](https://en.wikipedia.org/wiki/Fock_state) distribution of `ρ`. 
    
The `library` keyword argument specifies the plotting library to use, defaulting to [`Makie`](https://github.com/MakieOrg/Makie.jl). 

# Arguments
- `ρ::QuantumObject`: The quantum state for which to plot the Fock state distribution.
- `library::Union{Val,Symbol}`: The plotting library to use. Default is `Val(:Makie)`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. See the documentation for the specific plotting library for more information.

!!! note "Import library first"
    The plotting libraries must first be imported before using them with this function.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:Makie)` instead of `:Makie` as the plotting library. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
plot_fock_distribution(
    ρ::QuantumObject{SType};
    library::Union{Val,Symbol} = Val(:Makie),
    kwargs...,
) where {SType<:Union{Bra,Ket,Operator}} = plot_fock_distribution(makeVal(library), ρ; kwargs...)

plot_fock_distribution(::Val{T}, ρ::QuantumObject{SType}; kwargs...) where {T,SType<:Union{Bra,Ket,Operator}} =
    throw(ArgumentError("The specified plotting library $T is not available. Try running `using $T` first."))

@doc raw"""
    Bloch()

A structure representing a Bloch sphere visualization for quantum states."""
@kwdef mutable struct Bloch
    # Data storage
    """Points to plot on the Bloch sphere (3D coordinates)"""
    points::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}()
    """Vectors to plot on the Bloch sphere"""
    vectors::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    """Lines to draw on the sphere (points, style, properties)"""
    lines::Vector{Tuple{Vector{Vector{Float64}},String,Dict{Symbol,Any}}} =
        Vector{Tuple{Vector{Vector{Float64}},String,Dict{Symbol,Any}}}()
    """Arcs to draw on the sphere"""
    arcs::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}()

    # Style properties
    """Color of axis labels and text"""
    font_color::String = "#2E3440"
    """Font size for labels (default: 18)"""
    font_size::Int = 18
    """Transparency of wireframe"""
    frame_alpha::Float64 = 0.1
    """Color of wireframe"""
    frame_color::String = "#E5E9F0"
    """ Width of wireframe lines"""
    frame_width::Int = 1

    # Point properties
    """Default color cycle for points"""
    point_default_color::Vector{String} = ["blue", "red", "green", "orange", "cyan", "magenta", "yellow", "black"]
    """Colors for point markers"""
    point_color::Vector{String} = ["blue", "red", "green", "orange", "cyan", "magenta", "yellow", "black"]
    """Marker shapes (default: [:circle, :rect, :diamond, :utriangle])"""
    point_marker::Vector{Symbol} = [:circle, :rect, :diamond, :utriangle]
    """Marker sizes"""
    point_size::Vector{Int} = [40, 48, 50, 60]
    """Marker styles"""
    point_style::Vector{Symbol} = Symbol[]
    """Marker transparencies"""
    point_alpha::Vector{Float64} = Float64[]

    # Sphere properties
    """Transparency of Bloch sphere surface"""
    sphere_alpha::Float64 = 0.9
    """Color of Bloch sphere surface"""
    sphere_color::String = "#ECEFF4"

    # Layout properties
    """Figure size in pixels"""
    size::Tuple{Int,Int} = (700, 700)
    """Colors for vectors"""
    vector_color::Vector{String} = ["green", "blue", "orange", "red", "cyan", "magenta", "yellow", "black"]
    """Width of vectors (default: 2)"""
    vector_width::Int = 2
    """Azimuthal and elevation viewing angles in degrees (default: (-60, 30))"""
    view_angles::Tuple{Int,Int} = (-60, 30)

    # Label properties
    """Labels for x-axis (default: ["x", ""])"""
    xlabel::Vector{String} = ["x", ""]
    """Positions of x-axis labels"""
    xlpos::Vector{Float64} = [1.2, -1.2]
    """Labels for y-axis (default: ["y", ""])"""
    ylabel::Vector{String} = ["y", ""]
    """Positions of y-axis labels)"""
    ylpos::Vector{Float64} = [1.2, -1.2]
    """Labels for z-axis (default: ["|0⟩", "|1⟩"])"""
    zlabel::Vector{String} = ["|0⟩", "|1⟩"]
    """Positions of z-axis labels"""
    zlpos::Vector{Float64} = [1.2, -1.2]
end

@doc raw"""
    add_vectors!(b::Bloch, vec::Vector{<:Real})

Add a single normalized vector to the Bloch sphere visualization.

# Arguments
- `b::Bloch`: The Bloch sphere object to modify
- `vec::Vector{<:Real}`: A 3D vector to add (will be normalized)
- `vecs::Vector{<:Vector{<:Real}}}`: List of 3D vectors to add (each will be normalized)

# Example
```jldoctest
julia> b = Bloch();

julia> add_vectors!(b, [1, 0, 0])
1-element Vector{Vector{Float64}}:
 [1.0, 0.0, 0.0]
```

We can also add multiple normalized vectors to the Bloch sphere visualization.

```jldoctest
julia> b = Bloch();

julia> add_vectors!(b, [[1, 0, 0], [0, 1, 0]])
2-element Vector{Vector{Float64}}:
 [1.0, 0.0, 0.0]
 [0.0, 1.0, 0.0]
```
"""
function add_vectors!(b::Bloch, vec::Vector{<:Real})
    normalized_vec = normalize(convert(Vector{Float64}, vec))
    return push!(b.vectors, normalized_vec)
end
function add_vectors!(b::Bloch, vecs::Vector{<:Vector{<:Real}})
    return append!(b.vectors, [normalize(convert(Vector{Float64}, v)) for v in vecs])
end

@doc raw"""
    add_points!(b::Bloch, pnt::Vector{<:Real}; meth::Symbol = :s, color = nothing, alpha = 1.0)

Add a single point to the Bloch sphere visualization.

# Arguments
- b::Bloch: The Bloch sphere object to modify
- pnt::Vector{<:Real}: A 3D point to add
- meth::Symbol=:s: Display method (:s for single point, :m for multiple, :l for line)
- color: Color of the point (defaults to first default color if nothing)
- alpha=1.0: Transparency (1.0 = opaque, 0.0 = transparent)
"""
function add_points!(b::Bloch, pnt::Vector{<:Real}; meth::Symbol = :s, color = nothing, alpha = 1.0)
    points = reshape(convert(Vector{Float64}, pnt), 3, 1)
    return add_points!(b, points; meth = meth, color = color, alpha = alpha)
end

@doc raw"""
    add_points!(b::Bloch, pnts::Matrix{<:Real}; meth::Symbol = :s, color = nothing, alpha = 1.0)

Add multiple points to the Bloch sphere visualization.

# Arguments

- b::Bloch: The Bloch sphere object to modify
- pnts::Matrix{<:Real}: 3×N matrix of points (each column is a point)
- meth::Symbol=:s: Display method (:s for single point, :m for multiple, :l for line)
- color: Color of the points (defaults to first default color if nothing)
- alpha=1.0: Transparency (1.0 = opaque, 0.0 = transparent)
```
"""
function add_points!(b::Bloch, pnts::Matrix{<:Real}; meth::Symbol = :s, color = nothing, alpha = 1.0)
    if size(pnts, 1) != 3
        error("Points must be a 3×N matrix where columns are [x;y;z] points")
    end
    if meth ∉ [:s, :m, :l]
        error("meth must be :s, :m, or :l")
    end
    if meth == :s && size(pnts, 2) == 1
        pnts = hcat(pnts, pnts)
    end
    push!(b.points, convert(Matrix{Float64}, pnts))
    push!(b.point_style, meth)
    push!(b.point_alpha, alpha)
    return push!(b.point_color, color === nothing ? b.point_default_color[1] : color)
end

@doc raw"""
    add_line!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}; fmt = "k", kwargs...)

Add a line between two points on the Bloch sphere.

# Arguments
- b::Bloch: The Bloch sphere object to modify
- p1::Vector{<:Real}: First 3D point
- p2::Vector{<:Real}: Second 3D point
- fmt="k": Line format string (matplotlib style)
- kwargs...: Additional line properties
"""
function add_line!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}; fmt = "k", kwargs...)
    if length(p1) != 3 || length(p2) != 3
        error("Points must be 3D vectors")
    end
    x = [p1[2], p2[2]]
    y = [-p1[1], -p2[1]]
    z = [p1[3], p2[3]]
    return push!(b.lines, ([x, y, z], fmt, kwargs))
end

@doc raw"""
    add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real})

Add a circular arc through three points on the Bloch sphere.

# Arguments

- b::Bloch: The Bloch sphere object to modify
- p1::Vector{<:Real}: First 3D point
- p2::Vector{<:Real}: Second 3D point (middle point)
- p3::Vector{<:Real}: Third 3D point

# Examples

```jldoctest
julia> b = Bloch();

julia> add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
1-element Vector{Vector{Vector{Float64}}}:
 [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
```
"""
function add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real})
    return push!(b.arcs, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2)])
end
function add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real})
    return push!(b.arcs, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2), convert(Vector{Float64}, p3)])
end

@doc raw"""
    QuantumToolbox.add_states!(b::Bloch, states::QuantumObject...)

Add one or more quantum states to the Bloch sphere visualization by converting them into Bloch vectors.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify
- `states::QuantumObject...`: One or more quantum states (Ket, Bra, or Operator)

"""
function add_states! end

@doc raw"""
    clear!(b::Bloch)

Clear all graphical elements (points, vectors, lines, arcs) from the given Bloch sphere object `b`.

# Arguments

- `b::Bloch`  
  The Bloch sphere instance whose contents will be cleared.

# Returns

- The updated `Bloch` object `b` with all points, vectors, lines, and arcs removed.
"""
function clear!(b::Bloch)
    empty!(b.points)
    empty!(b.vectors)
    empty!(b.lines)
    empty!(b.arcs)
    return b
end

@doc raw"""
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
function render end

@doc raw"""
    plot_bloch(
        state::QuantumObject{<:Union{Ket,Bra,Operator}};
        library::Union{Symbol, Val} = :Makie,
        kwargs...
    )

Plot the state of a two-level quantum system on the Bloch sphere.

The `library` keyword argument specifies the plotting backend to use. The default is `:Makie`, which uses the [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) plotting library. This function internally dispatches to a type-stable version based on `Val(:Makie)` or other plotting backends.

# Arguments
- `state::QuantumObject`: The quantum state to be visualized. Can be a ket, bra, or operator.
- `library::Union{Symbol, Val}`: The plotting backend, either as a `Symbol` (e.g. `:Makie`) or a `Val` (e.g. `Val(:Makie)`). Default is `:Makie`.
- `kwargs...`: Additional keyword arguments passed to the specific plotting implementation.

!!! note "Import library first"
    The plotting backend library must be imported before use.

!!! warning "Beware of type-stability!"
    For improved performance and type-stability, prefer passing `Val(:Makie)` instead of `:Makie`. See [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) for details.
"""
function plot_bloch(state::QuantumObject{<:Union{Ket,Bra,Operator}}; library::Union{Symbol,Val} = :Makie, kwargs...)
    lib_val = library isa Symbol ? Val(library) : library
    return plot_bloch(lib_val, state; kwargs...)
end

function plot_bloch(::Val{T}, state::QuantumObject; kwargs...) where {T}
    return error("Unsupported backend: $T. Try :Makie or another supported library.")
end
