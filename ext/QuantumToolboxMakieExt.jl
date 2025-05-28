module QuantumToolboxMakieExt

using QuantumToolbox
using Makie:
    Axis, Axis3, Colorbar, Figure, GridLayout, heatmap!, surface!, barplot!, GridPosition, @L_str, Reverse, ylims!

include("bloch.jl")

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

"""
    _state_to_bloch(state::QuantumObject{<:Ket}) -> Vector{Float64}

Convert a quantum state (Ket) to its Bloch vector representation.

For a 2-level system (qubit), the Bloch vector components are calculated as:
r_i = ⟨ψ|σ_i|ψ⟩
where σ_i are the Pauli matrices.

For higher-dimensional systems, projects onto the generalized Bloch sphere.
"""
function _state_to_bloch(state::QuantumObject{<:Ket})
    if !isapprox(norm(state), 1.0, atol = 1e-6)
        @warn "State is not normalized. Normalizing before Bloch vector conversion."
        state = normalize(state)
    end
    N = length(state)
    if N == 2  # Qubit case
        ψ = state.data
        x = 2 * real(ψ[1] * conj(ψ[2]))
        y = 2 * imag(ψ[1] * conj(ψ[2]))
        z = abs2(ψ[1]) - abs2(ψ[2])
        return [x, y, z]
    else
        return _higher_dim_bloch_vector(state)
    end
end

"""
    _higher_dim_bloch_vector(state::QuantumObject{<:Ket}) -> Vector{Float64}

Compute the generalized Bloch vector for higher-dimensional systems using Gell-Mann basis.
"""
function _higher_dim_bloch_vector(state::QuantumObject{<:Ket})
    N = length(state)
    ψ = state.data
    # Number of generalized Bloch vector components: N²-1
    bloch_vec = zeros(Float64, N^2 - 1)
    # Symmetric (off-diagonal) components
    idx = 1
    for j in 1:N
        for k in (j+1):N
            bloch_vec[idx] = 2 * real(ψ[j] * conj(ψ[k]))
            bloch_vec[idx+1] = 2 * imag(ψ[j] * conj(ψ[k]))
            idx += 2
        end
    end
    # Diagonal components
    for l in 2:N
        for m in 1:(l-1)
            bloch_vec[idx] = sqrt(2/(l*(l-1))) * (abs2(ψ[m]) - (l-1)*abs2(ψ[l]))
            idx += 1
        end
    end
    return bloch_vec
end

"""
    plot_bloch(::Val{:Makie}, state::QuantumObject{<:Union{Ket,Bra}}; kwargs...)

Plot the state on a Bloch sphere using Makie.jl.

# Arguments

  - `state`: Quantum state to visualize (must be a Ket or Bra)
  - `show_axes`: Whether to show x/y/z axes (default: true)
  - `show_labels`: Whether to show axis labels (default: true)
  - `sphere_alpha`: Transparency of the sphere (default: 0.1)
  - `vector_color`: Color of the state vector (default: :red)
  - `kwargs...`: Additional arguments passed to the Bloch sphere renderer

# Returns

  - `fig`: The Makie Figure object
  - `ax`: The Axis3 object
  - `bloch`: The Bloch sphere object
"""
function QuantumToolbox.plot_bloch(::Val{:Makie}, state::QuantumObject{<:Union{Ket,Bra}}; kwargs...)
    state = isbra(state) ? dag(state) : state
    bloch_vec = _state_to_bloch(state)
    return _render_bloch_makie(bloch_vec; kwargs...)
end

"""
    _dm_to_bloch(ρ::QuantumObject{<:Operator}) -> Vector{Float64}

Convert a density matrix to its Bloch vector representation.
For qubits, uses Pauli matrices. For higher dimensions, uses generalized Gell-Mann basis.
"""
function _dm_to_bloch(ρ::QuantumObject{<:Operator})
    if !ishermitian(ρ)
        @warn "Density matrix is not Hermitian. Results may not be meaningful."
    end
    N = size(ρ, 1)

    if N == 2  # Qubit case
        σx = sigmax()
        σy = sigmay()
        σz = sigmaz()
        x = real(expect(σx, ρ))
        y = real(expect(σy, ρ))
        z = real(expect(σz, ρ))
        return [x, y, z]
    else
        return _higher_dim_bloch_vector(ρ)
    end
end

function QuantumToolbox.plot_bloch(::Val{:Makie}, ρ::QuantumObject{<:Operator}; kwargs...)
    bloch_vec = _dm_to_bloch(ρ)
    return _render_bloch_makie(bloch_vec; kwargs...)
end

function _render_bloch_makie(
    bloch_vec::Vector{Float64};
    location = nothing,
    show_axes = true,
    show_labels = true,
    sphere_alpha = 0.1,
    vector_color = :red,
    kwargs...,
)
    b = Bloch()
    b.sphere_alpha = sphere_alpha
    b.vector_color = [string(vector_color)]
    b.xlabel = show_labels ? ["x", "-x"] : ["", ""]
    b.ylabel = show_labels ? ["y", "-y"] : ["", ""]
    b.zlabel = show_labels ? ["|0⟩", "|1⟩"] : ["", ""]
    add_vectors!(b, bloch_vec)

    fig, location = _getFigAndLocation(location)
    fig, ax = render(b; location = location, kwargs...)
    return fig, ax
end

"""
    _higher_dim_bloch_vector(ρ::QuantumObject{<:Operator}) -> Vector{Float64}

Compute the generalized Bloch vector for higher-dimensional density matrices.
"""
function _higher_dim_bloch_vector(ρ::QuantumObject{<:Operator})
    N = size(ρ, 1)
    bloch_vec = zeros(Float64, N^2 - 1)
    # Symmetric (off-diagonal) components
    idx = 1
    for j in 1:N
        for k in (j+1):N
            bloch_vec[idx] = real(ρ[j, k] + ρ[k, j])
            bloch_vec[idx+1] = imag(ρ[j, k] - ρ[k, j])
            idx += 2
        end
    end
    # Diagonal components
    for l in 2:N
        for m in 1:(l-1)
            bloch_vec[idx] = sqrt(2/(l*(l-1))) * (real(ρ[m, m]) - (l-1)*real(ρ[l, l]))
            idx += 1
        end
    end
    return bloch_vec
end

end
