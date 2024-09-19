#=
Utilities:
    internal (or external) functions which will be used throughout the entire package
=#

export gaussian, n_th
export row_major_reshape, meshgrid

@doc raw"""
    row_major_reshape(Q::AbstractArray, shapes...)

Reshapes `Q` in the row-major order, as numpy.
"""
row_major_reshape(Q::AbstractArray{T}, shapes...) where {T} =
    PermutedDimsArray(reshape(Q, reverse(shapes)...), (length(shapes):-1:1))

@doc raw"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    return X, Y
end

@doc raw"""
    gaussian(x::Number, μ::Number, σ::Number)

Returns the gaussian function ``\exp \left[- \frac{(x - \mu)^2}{2 \sigma^2} \right]``,
where ``\mu`` and ``\sigma^2`` are the mean and the variance respectively.
"""
gaussian(x::Number, μ::Number, σ::Number) = exp(-(x - μ)^2 / (2 * σ^2))

@doc raw"""
    n_th(ω::Number, T::Real)

Gives the mean number of excitations in a mode with frequency ω at temperature T:
``n_{\rm th} (\omega, T) = \frac{1}{e^{\omega/T} - 1}``
"""
function n_th(ω::Real, T::Real)::Float64
    (T == 0 || ω == 0) && return 0.0
    abs(ω / T) > 50 && return 0.0
    return 1 / (exp(ω / T) - 1)
end

_get_dense_similar(A::AbstractArray, args...) = similar(A, args...)
_get_dense_similar(A::AbstractSparseMatrix, args...) = similar(nonzeros(A), args...)

_Ginibre_ensemble(n::Int, rank::Int = n) = randn(ComplexF64, n, rank) / sqrt(n)

makeVal(x::Val{T}) where {T} = x
makeVal(x) = Val(x)

getVal(x::Val{T}) where {T} = T

_get_size(A::AbstractMatrix) = size(A)
_get_size(A::AbstractVector) = (length(A), 1)

_non_static_array_warning(argname, arg::Tuple{}) =
    throw(ArgumentError("The argument $argname must be a Tuple or a StaticVector of non-zero length."))
_non_static_array_warning(argname, arg::Union{SVector{N,T},MVector{N,T},NTuple{N,T}}) where {N,T} = nothing
_non_static_array_warning(argname, arg::AbstractVector{T}) where {T} =
    @warn "The argument $argname should be a Tuple or a StaticVector for better performance. Try to use `$argname = $(Tuple(arg))` or `$argname = SVector(" *
          join(arg, ", ") *
          ")` instead of `$argname = $arg`." maxlog = 1

# convert tlist in time evolution
_convert_tlist(::Int32, tlist::AbstractVector) = _convert_tlist(Val(32), tlist)
_convert_tlist(::Float32, tlist::AbstractVector) = _convert_tlist(Val(32), tlist)
_convert_tlist(::ComplexF32, tlist::AbstractVector) = _convert_tlist(Val(32), tlist)
_convert_tlist(::Int64, tlist::AbstractVector) = _convert_tlist(Val(64), tlist)
_convert_tlist(::Float64, tlist::AbstractVector) = _convert_tlist(Val(64), tlist)
_convert_tlist(::ComplexF64, tlist::AbstractVector) = _convert_tlist(Val(64), tlist)
_convert_tlist(::Val{32}, tlist::AbstractVector) = convert(Vector{Float32}, tlist)
_convert_tlist(::Val{64}, tlist::AbstractVector) = convert(Vector{Float64}, tlist)