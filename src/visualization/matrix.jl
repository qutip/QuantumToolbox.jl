export matrix_histogram

@doc raw"""
    matrix_histogram(
        M::Union{QuantumObject,AbstractMatrix};
        library::Union{Val, Symbol} = Val(:Makie),
        method::Union{Symbol,Val} = Val(:real),
        kwargs...
    )

Plot a 3D histogram for the elements of matrix `M`.

The `library` keyword argument specifies the plotting library to use, defaulting to [`Makie`](https://github.com/MakieOrg/Makie.jl). 

# Arguments
- `M::Union{QuantumObject,AbstractMatrix}`: The [`QuantumObject`](@ref) or `AbstractMatrix` for which to be plotted. If it is a [`QuantumObject`](@ref), tt can be either a [`Operator`](@ref) or [`SuperOperator`](@ref).
- `library::Union{Val,Symbol}`: The plotting library to use. Default is `Val(:Makie)`.
- `method::Union{Symbol,Val}`: Method to use for plotting the matrix elements. Can be either `:real`, `:imag`, `:abs`, or `:angle`. Default is `Val(:real)`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function. See the documentation for the specific plotting library for more information.

!!! note "Import library first"
    The plotting libraries must first be imported before using them with this function.

!!! warning "Beware of type-stability!"
    If you want to keep type stability, it is recommended to use `Val(:Makie)` instead of `:Makie` as the plotting library. See [this link](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-value-type) and the [related Section](@ref doc:Type-Stability) about type stability for more details.
"""
matrix_histogram(
    M::Union{QuantumObject{QT},AbstractMatrix{MT}};
    library::Union{Val,Symbol} = Val(:Makie),
    method::Union{Symbol,Val} = Val(:real),
    kwargs...,
) where {QT<:Union{Operator,SuperOperator},MT<:Number} =
    matrix_histogram(makeVal(library), M; method = makeVal(method), kwargs...)

matrix_histogram(::Val{T}, M; kwargs...) where {T} =
    throw(ArgumentError("The specified plotting library $T is not available. Try running `using $T` first."))

# the following functions will be used in all plotting backends
_handle_matrix_plot_data(M::QuantumObject{T}, method::Val) where {T<:Union{Operator,SuperOperator}} =
    _handle_matrix_plot_data(M.data, method)
_handle_matrix_plot_data(M::AbstractMatrix{T}, ::Val{:real}) where {T<:Number} = real(transpose(M))
_handle_matrix_plot_data(M::AbstractMatrix{T}, ::Val{:imag}) where {T<:Number} = imag(transpose(M))
_handle_matrix_plot_data(M::AbstractMatrix{T}, ::Val{:abs}) where {T<:Number} = abs.(transpose(M))
_handle_matrix_plot_data(M::AbstractMatrix{T}, ::Val{:angle}) where {T<:Number} = angle.(transpose(M))
_handle_matrix_plot_data(::AbstractMatrix{T}, method::Val) where {T<:Number} =
    throw(ArgumentError("Invalid keyword argument method = $(method), should be either: :real, :imag, :abs, or :angle"))

# for y-axis ticks
_gen_default_ket_labels(::QuantumObject{SuperOperator}, ydata) = map(y -> L"|%$(y)\rangle\!\rangle", ydata)
_gen_default_ket_labels(::Union{QuantumObject{Operator},AbstractMatrix{T}}, ydata) where {T<:Number} =
    map(y -> L"|%$(y)\rangle", ydata)

# for x-axis ticks
_gen_default_bra_labels(::QuantumObject{SuperOperator}, xdata) = map(x -> L"\langle\!\langle%$(x)|", xdata)
_gen_default_bra_labels(::Union{QuantumObject{Operator},AbstractMatrix{T}}, xdata) where {T<:Number} =
    map(x -> L"\langle%$(x)|", xdata)