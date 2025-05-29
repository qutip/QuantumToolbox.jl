export plot_wigner,
    plot_fock_distribution, plot_bloch, Bloch, add_points!, add_vectors!, add_line!, add_arc!, clear!, render

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

mutable struct Bloch
    points::Vector{Vector{Float64}}
    vectors::Vector{Vector{Float64}}
    lines::Vector{Vector{Vector{Float64}}}
    arcs::Vector{Vector{Vector{Float64}}}
    font_color::String
    font_size::Int
    frame_alpha::Float64
    frame_color::String
    frame_width::Int
    point_color::Vector{String}
    point_marker::Vector{Symbol}
    point_size::Vector{Int}
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
            "black",           # Font color
            18,                # Slightly smaller font size
            0.0,               # Transparent frame (minimal background)
            "white",           # Background color (invisible if alpha = 0)
            0,                 # No visible frame grid
            ["blue", "red"],   # Fewer point colors
            [:circle, :rect],  # Simplified marker shapes
            [6, 6],            # Uniform point sizes
            0.8,               # High opacity sphere
            "#DDDDFF",         # Subtle light-blue sphere
            [800, 800],        # Slightly larger canvas
            ["blue", "red"],   # Fewer vector colors
            2,                 # Thinner vector lines
            [-60, 30],         # Standard Bloch viewing angle
            ["x", "-x"],
            [1.05, -1.05],     # Slightly closer label position
            ["y", "-y"],
            [1.05, -1.05],
            ["|0⟩", "|1⟩"],
            [1.1, -1.1],
            nothing,
            nothing,
        )
    end
end

function add_points!(b::Bloch, pnt::Vector{<:Real}) 
    push!(b.points, convert(Vector{Float64}, pnt))
end

function add_points!(b::Bloch, pnts::Vector{<:Vector{<:Real}}) 
    append!(b.points, [convert(Vector{Float64}, p) for p in pnts])
end

function add_vectors!(b::Bloch, vec::Vector{<:Real})
    normalized_vec = normalize(convert(Vector{Float64}, vec))
    push!(b.vectors, normalized_vec)
end

function add_vectors!(b::Bloch, vecs::Vector{<:Vector{<:Real}})
    append!(b.vectors, [normalize(convert(Vector{Float64}, v)) for v in vecs])
end

function add_line!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real})
    push!(b.lines, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2)])
end

function add_arc!(b::Bloch, p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real})
    push!(b.arcs, [convert(Vector{Float64}, p1), convert(Vector{Float64}, p2), convert(Vector{Float64}, p3)])
end

function clear!(b::Bloch)
    empty!(b.points)
    empty!(b.vectors)
    empty!(b.lines)
    empty!(b.arcs)
    return b
end

function render end

function plot_bloch(state::QuantumObject{<:Union{Ket,Bra,Operator}}; library::Union{Symbol,Val}=:Makie, kwargs...)
    lib_val = library isa Symbol ? Val(library) : library
    return plot_bloch(lib_val, state; kwargs...)
end

function plot_bloch(::Val{T}, state::QuantumObject; kwargs...) where {T}
    error("Unsupported backend: $T. Try :Makie or another supported library.")
end