export plot_wigner
export plot_fock_distribution
export plot_bloch, Bloch, render, add_points!, add_vectors!, add_line!, add_arc!, clear!, add_states!

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
    Bloch(kwargs...)

A structure representing a Bloch sphere visualization for quantum states. Available keyword arguments are listed in the following fields.

# Fields:

## Data storage
- `points::Vector{Matrix{Float64}}`: Points to plot on the Bloch sphere (3D coordinates)
- `vectors::Vector{Vector{Float64}}}`: Vectors to plot on the Bloch sphere
- `lines::Vector{Tuple{Vector{Vector{Float64}},String,Dict{Any,Any}}}`: Lines to draw on the sphere (points, style, properties)
- `arcs::Vector{Vector{Vector{Float64}}}}`: Arcs to draw on the sphere

## Style properties

- `font_color::String`: Color of axis labels and text
- `font_size::Int`: Font size for labels. Default: `20`
- `frame_alpha::Float64`: Transparency of the wire frame
- `frame_color::String`: Color of the wire frame

## Point properties

- `point_default_color::Vector{String}}`: Default color cycle for points
- `point_color::Vector{String}}`: List of colors for Bloch point markers to cycle through
- `point_marker::Vector{Symbol}}`: List of point marker shapes to cycle through. Default: `[:circle, :rect, :diamond, :utriangle]`
- `point_size::Vector{Int}}`: List of point marker sizes (not all markers look the same size when plotted)
- `point_style::Vector{Symbol}}`: List of marker styles
- `point_alpha::Vector{Float64}}`: List of marker transparencies

## Sphere properties

- `sphere_color::String`: Color of Bloch sphere surface
- `sphere_alpha::Float64`: Transparency of sphere surface. Default: `0.2`

## Vector properties

- `vector_color`::Vector{String}: Colors for vectors
- `vector_width`::Float64: Width of vectors
- `vector_arrowsize`::Vector{Float64}: Arrow size parameters as [head length, head width, stem width]

## Layout properties

- `view::Vector{Int}`: Azimuthal and elevation viewing angles in degrees. Default: `[30, 30]`

## Label properties

- `xlabel::Vector{AbstractString}`: Labels for x-axis. Default: `[L"x", ""]`
- `xlpos::Vector{Float64}`: Positions of x-axis labels. Default: `[1.2, -1.2]`
- `ylabel::Vector{AbstractString}`: Labels for y-axis. Default: `[L"y", ""]`
- `ylpos::Vector{Float64}`: Positions of y-axis labels. Default: `[1.2, -1.2]`
- `zlabel::Vector{AbstractString}`: Labels for z-axis. Default: `[L"|0\rangle", L"|1\rangle"]`
- `zlpos::Vector{Float64}`: Positions of z-axis labels. Default: `[1.2, -1.2]`
"""
@kwdef mutable struct Bloch
    points::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}()
    vectors::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    lines::Vector{Tuple{Vector{Vector{Float64}},String}} = Vector{Tuple{Vector{Vector{Float64}},String}}()
    arcs::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}()
    font_color::String = "black"
    font_size::Int = 20
    frame_alpha::Float64 = 0.1
    frame_color::String = "gray"
    point_default_color::Vector{String} = ["blue", "red", "green", "#CC6600"]
    point_color::Vector{Union{Nothing,String}} = Union{Nothing,String}[]
    point_marker::Vector{Symbol} = [:circle, :rect, :diamond, :utriangle]
    point_size::Vector{Float64} = [5.5, 6.2, 6.5, 7.5]
    point_style::Vector{Symbol} = Symbol[]
    point_alpha::Vector{Float64} = Float64[]
    sphere_alpha::Float64 = 0.2
    sphere_color::String = "#FFDDDD"
    vector_color::Vector{String} = ["green", "#CC6600", "blue", "red"]
    vector_width::Float64 = 0.025
    vector_arrowsize::Vector{Float64} = [0.07, 0.08, 0.08]
    view::Vector{Int} = [30, 30]
    xlabel::Vector{AbstractString} = [L"x", ""]
    xlpos::Vector{Float64} = [1.2, -1.2]
    ylabel::Vector{AbstractString} = [L"y", ""]
    ylpos::Vector{Float64} = [1.2, -1.2]
    zlabel::Vector{AbstractString} = [L"|0\rangle", L"|1\rangle"]
    zlpos::Vector{Float64} = [1.2, -1.2]
end

const BLOCH_DATA_FIELDS = (:points, :vectors, :lines, :arcs)
function Base.show(io::IO, b::Bloch)
    # To align the output and make it easier to read
    # we use rpad `17` and `19` for Bloch sphere data and properties, respectively
    # 17 is the length of string: `Number of vectors`
    # 19 is the length of string: `point_default_color`
    println(io, "Bloch Sphere\n")
    println(io, "data:")
    println(io, "-----")
    map(n -> println(io, rpad("Number of $n", 17, " "), " = ", length(getfield(b, n))), BLOCH_DATA_FIELDS)
    println(io, "")
    println(io, "properties:")
    println(io, "-----------")
    map(n -> (n ∉ BLOCH_DATA_FIELDS) && (println(io, rpad("$n", 19, " "), " = ", getfield(b, n))), fieldnames(Bloch))
    return nothing
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
add_vectors!(b::Bloch, vec::Vector{<:Real}) = push!(b.vectors, convert(Vector{Float64}, vec))
add_vectors!(b::Bloch, vecs::Vector{<:Vector{<:Real}}) = append!(b.vectors, [convert(Vector{Float64}, v) for v in vecs])

@doc raw"""
    add_points!(b::Bloch, pnt::Vector{<:Real}; meth::Symbol = :s, color = "blue", alpha = 1.0)

Add a single point to the Bloch sphere visualization.

# Arguments
- `b::Bloch`: The Bloch sphere object to modify
- `pnt::Vector{Float64}`: A 3D point to add
- `meth::Symbol=:s`: Display method (`:s` for single point, `:m` for multiple, `:l` for line)
- `color`: Color of the point (defaults to first default color if nothing)
- `alpha=1.0`: Transparency (`1.0` means opaque and `0.0` means transparent)
"""
function add_points!(b::Bloch, pnt::Vector{Float64}; meth::Symbol = :s, color = nothing, alpha = 1.0)
    return add_points!(b, reshape(pnt, 3, 1); meth, color, alpha)
end
function add_points!(b::Bloch, pnts::Vector{Vector{Float64}}; meth::Symbol = :s, color = nothing, alpha = 1.0)
    return add_points!(b, Matrix(hcat(pnts...)'); meth, color, alpha)
end

@doc raw"""
    add_points!(b::Bloch, pnts::Matrix{Float64}; meth::Symbol = :s, color = nothing, alpha = 1.0)

Add multiple points to the Bloch sphere visualization.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify
- `pnts::Matrix{Float64}`: `3×N` matrix of points (each column is a point)
- `meth::Symbol=:s`: Display method (`:s` for single point, `:m` for multiple, `:l` for line)
- `color`: Color of the points (defaults to first default color if nothing)
- `alpha=1.0`: Transparency (`1.0` means opaque and `0.0` means transparent)
```
"""
function add_points!(
    b::Bloch,
    pnts::Matrix{<:Real};
    meth::Symbol = :s,
    color::Union{Nothing,String} = nothing,
    alpha::Float64 = 1.0,
)
    (size(pnts, 1) == 3) || throw(ArgumentError("Points must be a 3×N matrix where each column is [x; y; z]"))
    (meth in (:s, :m, :l)) || throw(ArgumentError("`meth` must be :s, :m, or :l"))

    push!(b.points, convert(Matrix{Float64}, pnts))
    push!(b.point_style, meth)
    push!(b.point_alpha, alpha)
    push!(b.point_color, color)
    return nothing
end

@doc raw"""
    add_line!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}; fmt = "k", kwargs...)

Add a line between two points on the Bloch sphere.

# Arguments
- `b::Bloch`: The Bloch sphere object to modify
- `p1::Vector{<:Real}`: First 3D point
- `p2::Vector{<:Real}`: Second 3D point
- `fmt="k"`: Line format string (matplotlib style)
"""
function add_line!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}; fmt = "k")
    (length(p1) != 3 || length(p2) != 3) && throw(ArgumentError("Points must be 3D vectors"))
    x = [p1[1], p2[1]]
    y = [p1[2], p2[2]]
    z = [p1[3], p2[3]]
    push!(b.lines, ([x, y, z], fmt))
    return b
end

@doc raw"""
    add_line!(
        b::Bloch,
        start_point::QuantumObject,
        end_point::QuantumObject;
        fmt = "k"
    )

Add a line between two quantum states on the Bloch sphere visualization.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify.
- `start_point::QuantumObject`: The starting quantum state. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `end_point::QuantumObject`: The ending quantum state. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `fmt::String="k"`: (optional) A format string specifying the line style and color (default is black `"k"`).

# Description

This function converts the given quantum states into their Bloch vector representations and adds a line between these two points on the Bloch sphere visualization. 

# Example

```julia
b = Bloch()
ψ₁ = normalize(basis(2, 0) + basis(2, 1))
ψ₂ = normalize(basis(2, 0) - im * basis(2, 1))
add_line!(b, ψ₁, ψ₂; fmt = "r--")
```
"""
function add_line!(
    b::Bloch,
    start_point::QuantumObject{OpType1},
    end_point::QuantumObject{OpType2};
    fmt = "k",
) where {OpType1<:Union{Ket,Bra,Operator},OpType2<:Union{Ket,Bra,Operator}}
    coords1 = _state_to_bloch(start_point)
    coords2 = _state_to_bloch(end_point)
    return add_line!(b, coords1, coords2; fmt = fmt)
end

@doc raw"""
    add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real})

Add a circular arc through three points on the Bloch sphere.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify
- `p1::Vector{<:Real}`: Starting 3D point
- `p2::Vector{<:Real}`: [Optional] Middle 3D point
- `p3::Vector{<:Real}`: Ending 3D point

# Examples

```jldoctest
julia> b = Bloch();

julia> add_arc!(b, [1, 0, 0], [0, 1, 0], [0, 0, 1])
1-element Vector{Vector{Vector{Float64}}}:
 [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
```
"""
function add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real})
    (length(p1) != 3 || length(p2) != 3) && throw(ArgumentError("Points must be 3D vectors"))
    return push!(b.arcs, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2)])
end
function add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real})
    (length(p1) != 3 || length(p2) != 3 || length(p3) != 3) && throw(ArgumentError("Points must be 3D vectors"))
    return push!(b.arcs, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2), convert(Vector{Float64}, p3)])
end

@doc raw"""
    add_arc!(
        b::Bloch,
        start_point::QuantumObject,
        middle_point::QuantumObject,
        end_point::QuantumObject
    )

Add a circular arc through three points on the Bloch sphere.

# Arguments

- `b::Bloch`: The Bloch sphere object to modify.
- `start_point::QuantumObject`: The starting quantum state. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `middle_point::QuantumObject`: [Optional] The middle quantum state. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `end_point::QuantumObject`: The ending quantum state. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).

# Description

This function converts the given quantum states into their Bloch vector representations and adds a arc between these two (or three) points on the Bloch sphere visualization. 
"""
function add_arc!(
    b::Bloch,
    start_point::QuantumObject{OpType1},
    end_point::QuantumObject{OpType2},
) where {OpType1<:Union{Ket,Bra,Operator},OpType2<:Union{Ket,Bra,Operator}}
    coords1 = _state_to_bloch(start_point)
    coords2 = _state_to_bloch(end_point)
    return add_arc!(b, coords1, coords2)
end
function add_arc!(
    b::Bloch,
    start_point::QuantumObject{OpType1},
    middle_point::QuantumObject{OpType2},
    end_point::QuantumObject{OpType3},
) where {OpType1<:Union{Ket,Bra,Operator},OpType2<:Union{Ket,Bra,Operator},OpType3<:Union{Ket,Bra,Operator}}
    coords1 = _state_to_bloch(start_point)
    coords2 = _state_to_bloch(middle_point)
    coords3 = _state_to_bloch(end_point)
    return add_arc!(b, coords1, coords2, coords3)
end

@doc raw"""
    add_states!(b::Bloch, states::Vector{QuantumObject})

Add one or more quantum states to the Bloch sphere visualization by converting them into Bloch vectors.

# Arguments
- `b::Bloch`: The Bloch sphere object to modify
- `states::Vector{QuantumObject}`: One or more quantum states ([`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref))

# Example

```julia
x = basis(2, 0) + basis(2, 1);
y = basis(2, 0) + im * basis(2, 1);
z = basis(2, 0);
b = Bloch();
add_states!(b, [x, y, z])
```
"""
function add_states!(b::Bloch, states::Vector{<:QuantumObject})
    vecs = map(state -> _state_to_bloch(state), states)
    append!(b.vectors, vecs)
    return b.vectors
end

function add_states!(b::Bloch, state::QuantumObject)
    push!(b.vectors, _state_to_bloch(state))
    return b.vectors
end

_state_to_bloch(state::QuantumObject{Ket}) = _ket_to_bloch(state)
_state_to_bloch(state::QuantumObject{Bra}) = _ket_to_bloch(state')
_state_to_bloch(state::QuantumObject{Operator}) = _dm_to_bloch(state)

raw"""
    _ket_to_bloch(state::QuantumObject{Ket}) -> Vector{Float64}

Convert a pure qubit state (`Ket`) to its Bloch vector representation.

If the state is not normalized, it is automatically normalized before conversion.

# Arguments
- `state`: A `Ket` representing a pure quantum state.

# Returns
A 3-element `Vector{Float64}` representing the Bloch vector `[x, y, z]`.

# Throws
- `ArgumentError` if the state dimension is not 2.
"""
function _ket_to_bloch(state::QuantumObject{Ket})
    (size(state) == (2,)) ||
        throw(ArgumentError("Bloch sphere visualization is only supported for qubit states (2-level systems)"))

    state_norm = norm(state)
    if !isapprox(state_norm, 1.0, atol = 1e-6)
        @warn "State is not normalized. Normalizing before Bloch vector conversion."
        ψ = state.data / state_norm
    else
        ψ = state.data
    end

    c = conj(ψ[1]) * ψ[2]
    x = 2 * real(c)
    y = 2 * imag(c)
    z = abs2(ψ[1]) - abs2(ψ[2])
    return [x, y, z]
end

raw"""
    _dm_to_bloch(ρ::QuantumObject{Operator}) -> Vector{Float64}

Convert a qubit density matrix (`Operator`) to its Bloch vector representation.

This function assumes the input is Hermitian. If the density matrix is not Hermitian, a warning is issued.

# Arguments
- `ρ`: A density matrix (`Operator`) representing a mixed or pure quantum state.

# Returns
A 3-element `Vector{Float64}` representing the Bloch vector `[x, y, z]`.

# Throws
- `ArgumentError` if the matrix dimension is not 2.
"""
function _dm_to_bloch(ρ::QuantumObject{Operator})
    (size(ρ) == (2, 2)) ||
        throw(ArgumentError("Bloch sphere visualization is only supported for qubit states (2-level systems)"))

    ishermitian(ρ) || (@warn "Density matrix is not Hermitian. Results may not be meaningful.")

    state_norm = norm(ρ)
    if !isapprox(state_norm, 1.0, atol = 1e-6)
        @warn "State is not normalized. Normalizing before Bloch vector conversion."
        ρ2 = ρ / state_norm
    else
        ρ2 = ρ
    end
    x = real(ρ2[1, 2] + ρ2[2, 1])
    y = imag(ρ2[2, 1] - ρ2[1, 2])
    z = real(ρ2[1, 1] - ρ2[2, 2])
    return [x, y, z]
end

@doc raw"""
    clear!(b::Bloch)

Clear all graphical elements (points, vectors, lines, arcs) from the given [`Bloch`](@ref) sphere object `b`.

# Arguments

- `b::Bloch`: The Bloch sphere instance whose contents will be cleared.

# Returns

- The updated `Bloch` object `b` with all points, vectors, lines, and arcs removed.
"""
function clear!(b::Bloch)
    empty!(b.points)
    empty!(b.point_color)
    empty!(b.point_style)
    empty!(b.point_alpha)
    empty!(b.vectors)
    empty!(b.lines)
    empty!(b.arcs)
    return b
end

@doc raw"""
    render(b::Bloch; location=nothing)

Render the Bloch sphere visualization from the given [`Bloch`](@ref) object `b`.

# Arguments

- `b::Bloch`: The Bloch sphere object containing states, vectors, and settings to visualize.
- `location::Union{GridPosition,Nothing}`: The location of the plot in the layout. If `nothing`, the plot is created in a new figure. Default is `nothing`.

# Returns

- A tuple `(fig, lscene)` where `fig` is the figure object and `lscene` is the `LScene` object used for plotting. These can be further manipulated or saved by the user.
"""
function render end

@doc raw"""
    plot_bloch(
        state::QuantumObject;
        library::Union{Symbol, Val} = :Makie,
        kwargs...
    )

Plot the state of a two-level quantum system on the Bloch sphere.

The `library` keyword argument specifies the plotting backend to use. The default is `:Makie`, which uses the [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) plotting library. This function internally dispatches to a type-stable version based on `Val(:Makie)` or other plotting backends.

# Arguments
- `state::QuantumObject`: The quantum state to be visualized. Can be a [`Ket`](@ref), [`Bra`](@ref), or [`Operator`](@ref).
- `library::Union{Val,Symbol}`: The plotting library to use. Default is `Val(:Makie)`.
- `kwargs...`: Additional keyword arguments passed to the specific plotting implementation.

!!! note "Import library first"
    The plotting backend library must be imported before use.

!!! warning "Beware of type-stability!"
    For improved performance and type-stability, prefer passing `Val(:Makie)` instead of `:Makie`. See [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) for details.
"""
plot_bloch(
    state::QuantumObject{OpType};
    library::Union{Symbol,Val} = Val(:Makie),
    kwargs...,
) where {OpType<:Union{Ket,Bra,Operator}} = plot_bloch(makeVal(library), state; kwargs...)
plot_bloch(::Val{T}, state::QuantumObject{OpType}; kwargs...) where {T,OpType<:Union{Ket,Bra,Operator}} =
    throw(ArgumentError("The specified plotting library $T is not available. Try running `using $T` first."))
