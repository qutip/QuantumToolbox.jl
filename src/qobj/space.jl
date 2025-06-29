#=
This file defines the Hilbert space structure.
=#

export AbstractSpace, Space

abstract type AbstractSpace end

@doc raw"""
    struct Space <: AbstractSpace
        size::Int
    end

A structure that describes a single Hilbert space with size = `size`.
"""
struct Space <: AbstractSpace
    size::Int

    function Space(size::Int)
        (size < 1) && throw(DomainError(size, "The size of Space must be positive integer (≥ 1)."))
        return new(size)
    end
end

dimensions_to_dims(s::Space) = SVector{1,Int}(s.size)

# this creates a list of Space(1), it is used to generate `from` for Ket, and `to` for Bra)
space_one_list(N::Int) = ntuple(i -> Space(1), Val(N))
