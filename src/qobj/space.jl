export AbstractSpace, Field, Space

abstract type AbstractSpace end

# this replaces Space(1), so that we don't need to store the value `1`
struct Field <: AbstractSpace end
Base.getproperty(s::Field, key::Symbol) = getproperty(s, Val{key}())
Base.getproperty(s::Field, ::Val{:size}) = 1

# this creates a list of Field, it is used to generate `from` for Ket, and `to` for Bra)
Field_list(N::Int) = SVector{N,AbstractSpace}(ntuple(i -> Field(), Val(N)))

struct Space <: AbstractSpace
    size::Int

    function Space(size::Int)
        # (put this comment here to avoid bad JuliaFormatter syntax)
        if size > 1
            return new(size)
        elseif size == 1
            return Field()
        else
            throw(DomainError(size, "The size of Space must be positive integer (≥ 1)."))
        end
    end
end

dims_to_list(s::SType) where {SType<:Union{Field,Space}} = SVector{1,Int}(s.size)

# TODO: introduce energy restricted space
#=
struct EnrSpace{N} <: AbstractSpace
    size::Int
    dims::SVector{N,Int}
    n_excitations::Int
    state2idx
    idx2state
end

dims_to_list(s::EnrSpace) = s.dims
=#
