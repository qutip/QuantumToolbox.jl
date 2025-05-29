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

A structure representing a Bloch sphere visualization for quantum states.

# Fields
- `points::Vector{Matrix{Float64}}`: Points to plot on the Bloch sphere (3D coordinates)
- `vectors::Vector{Vector{Float64}}`: Vectors to plot on the Bloch sphere
- `lines::Vector{Tuple{Vector{Vector{Float64}},String,Dict{Any,Any}}}`: Lines to draw on the sphere (points, style, properties)
- `arcs::Vector{Vector{Vector{Float64}}}`: Arcs to draw on the sphere
- `font_color::String`: Color of axis labels and text (default: "#2E3440")
- `font_size::Int`: Font size for labels (default: 18)
- `frame_alpha::Float64`: Transparency of wireframe (default: 0.1)
- `frame_color::String`: Color of wireframe (default: "#E5E9F0")
- `frame_width::Int`: Width of wireframe lines (default: 1)
- `point_default_color::Vector{String}`: Default color cycle for points (default: ["blue", "red", "green", "orange"])
- `point_color::Vector{Any}`: Colors for point markers (default: ["blue", "red", "green", "orange"])
- `point_marker::Vector{Symbol}`: Marker shapes (default: [:circle, :rect, :diamond, :utriangle])
- `point_size::Vector{Int}`: Marker sizes (default: [40, 48, 50, 60])
- `point_style::Vector{Symbol}`: Marker styles (default: empty)
- `point_alpha::Vector{Float64}`: Marker transparencies (default: empty)
- `sphere_alpha::Float64`: Transparency of Bloch sphere surface (default: 0.8)
- `sphere_color::String`: Color of Bloch sphere surface (default: "#ECEFF4")
- `size::Vector{Int}`: Figure size in pixels (default: [800, 800])
- `vector_color::Vector{String}`: Colors for vectors (default: ["green", "blue", "orange"])
- `vector_width::Int`: Width of vectors (default: 2)
- `view_angles::Vector{Int}`: Azimuthal and elevation viewing angles in degrees (default: [-60, 30])
- `xlabel::Vector{String}`: Labels for x-axis (default: ["x", ""])
- `xlpos::Vector{Float64}`: Positions of x-axis labels (default: [1.2, -1.2])
- `ylabel::Vector{String}`: Labels for y-axis (default: ["y", ""])
- `ylpos::Vector{Float64}`: Positions of y-axis labels (default: [1.2, -1.2])
- `zlabel::Vector{String}`: Labels for z-axis (default: ["|0⟩", "|1⟩"])
- `zlpos::Vector{Float64}`: Positions of z-axis labels (default: [1.2, -1.2])
- `fig::Any`: Figure object (default: nothing)
- `ax::Any`: Axes object (default: nothing)
"""
mutable struct Bloch
    points::Vector{Matrix{Float64}}
    vectors::Vector{Vector{Float64}}
    lines::Vector{Tuple{Vector{Vector{Float64}},String,Dict{Any,Any}}}
    arcs::Vector{Vector{Vector{Float64}}}
    font_color::String
    font_size::Int
    frame_alpha::Float64
    frame_color::String
    frame_width::Int
    point_default_color::Vector{String}
    point_color::Vector{Any}
    point_marker::Vector{Symbol}
    point_size::Vector{Int}
    point_style::Vector{Symbol}
    point_alpha::Vector{Float64}
    sphere_alpha::Float64
    sphere_color::String
    size::Vector{Int}
    vector_color::Vector{String}
    vector_width::Int
    view_angles::Vector{Int}
    xlabel::Vector{String}
    xlpos::Vector{Float64}
    ylabel::Vector{String}
    ylpos::Vector{Float64}
    zlabel::Vector{String}
    zlpos::Vector{Float64}
    fig::Any
    ax::Any

    function Bloch()
        return new(
            [],
            [],
            [],
            [],
            "#2E3440",
            18,
            0.1,
            "#E5E9F0",
            1,
            ["blue", "red", "green", "orange"],
            ["blue", "red", "green", "orange"],
            [:circle, :rect, :diamond, :utriangle],
            [40, 48, 50, 60],
            [],
            [],
            0.8,
            "#ECEFF4",
            [800, 800],
            ["green", "blue", "orange"],
            2,
            [-60, 30],
            ["x", ""],
            [1.2, -1.2],
            ["y", ""],
            [1.2, -1.2],
            ["|0⟩", "|1⟩"],
            [1.2, -1.2],
            nothing,
            nothing,
        )
    end
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

# Examples

```jldoctest
julia> b = Bloch();

julia> add_line!(b, [1, 0, 0], [0, 1, 0])
1-element Vector{Tuple{Vector{Vector{Float64}}, String, Dict{Any, Any}}}:
 ([[0.0, 1.0], [-1.0, 0.0], [0.0, 0.0]], "k", Dict())
```
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
