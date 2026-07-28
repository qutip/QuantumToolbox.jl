export meshgrid

@doc raw"""
    meshgrid(x::AbstractVector, y::AbstractVector)

Equivalent to [numpy meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html).
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    X = reshape(repeat(x, inner = length(y)), length(y), length(x))
    Y = repeat(y, outer = (1, length(x)))
    return X, Y
end
