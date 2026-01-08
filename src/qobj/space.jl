#=
This file defines the Hilbert space structure.
=#

export AbstractSpace, HilbertSpace

@doc raw"""
    abstract type AbstractSpace

Abstract type for all Hilbert and Liouville space structures.
"""
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
        (size < 1) && throw(DomainError(size, "The size of `HilbertSpace` must be positive integer (≥ 1)."))
        return new(size)
    end
end

Base.length(s::HilbertSpace) = 1

Base.length(s::HilbertSpace) = 1

dimensions_to_dims(s::HilbertSpace) = SVector{1, Int}(s.size)

get_hilbert_size(s::HilbertSpace) = s.size
get_liouville_size(s::HilbertSpace) = s.size^2

# this creates a list of HilbertSpace(1), it is used to generate `from` for Ket, and `to` for Bra
hilbertspace_one_list(dimensions::NTuple{N,AbstractSpace}) where {N} =
    ntuple(i -> HilbertSpace(1), Val(sum(length, dimensions)))
