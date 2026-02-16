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
        fock_numbers::Union{Nothing, AbstractVector} = nothing,
        unit_y_range::Bool = true,
        location::Union{GridPosition, Nothing} = nothing,
        kwargs...,
    ) where {SType <: Union{Bra, Ket, Operator}}
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
        fock_numbers::Union{Nothing, AbstractVector} = nothing,
        unit_y_range::Bool = true,
        location::Union{GridPosition, Nothing} = nothing,
        kwargs...,
    ) where {SType <: Union{Bra, Ket, Operator}}
    ρ = ket2dm(ρ)
    D = get_size(ρ.dimensions)[1]
    isapprox(tr(ρ), 1, atol = 1.0e-4) || (@warn "The input ρ should be normalized.")

    # handle x ticks
    xvec = 0:(D - 1)
    fock_numbers = isnothing(fock_numbers) ? string.(collect(xvec)) : fock_numbers
    length(fock_numbers) == D ||
        throw(ArgumentError("Length of fock_numbers ($(length(fock_numbers))) does not match the total dimension: $D"))

    # handle limits
    # -0.5 originates from the half width of bar
    limits = unit_y_range ? (-0.5, D - 0.5, 0, 1) : (-0.5, D - 0.5, nothing, nothing)

    fig, location = _getFigAndLocation(location)
    lyt = GridLayout(location)
    ax = Axis(
        lyt[1, 1];
        xticks = (xvec, fock_numbers),
        xlabel = "Fock number",
        ylabel = "Occupation probability",
        limits = limits,
    )

    bp = barplot!(ax, xvec, real(diag(ρ)); kwargs...)

    return fig, ax, bp
end
