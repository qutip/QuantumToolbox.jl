#=
This file defines the Hilbert space structure.
=#

export AbstractSpace, Space

abstract type AbstractSpace end

@doc raw"""
    struct HilbertSpace <: AbstractSpace
        size::Int
    end

A structure that describes a single Hilbert space with size = `size`.
"""
struct HilbertSpace <: AbstractSpace
    size::Int

    function HilbertSpace(size::Int)
        (size < 1) && throw(DomainError(size, "The size of Space must be positive integer (≥ 1)."))
        return new(size)
    end
end

dimensions_to_dims(s::HilbertSpace) = SVector{1,Int}(s.size)
hilbert_dimensions_to_size(s::HilbertSpace) = s.size
liouvillian_dimensions_to_size(s::HilbertSpace) = s.size^2
