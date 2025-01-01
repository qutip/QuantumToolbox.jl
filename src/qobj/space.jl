export AbstractSpace, Field, Space

abstract type AbstractSpace end

# this replaces Space(1), so that we don't need to store the value `1`
struct Field <: AbstractSpace end
Base.getproperty(s::Field, key::Symbol) = getproperty(s, Val{key}())
Base.getproperty(s::Field, ::Val{:size}) = 1
Base.string(s::Field) = "1" # this is only used when printing AbstractDimensions

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
Base.string(s::Space) = string(s.size) # this is only used when printing AbstractDimensions

# for `prod(::Dimensions)`
Base.:(*)(i::Int, s::AbstractSpace) = i * s.size
Base.:(*)(s1::AbstractSpace, s2::AbstractSpace) = s1.size * s2.size
Base.prod(spaces::SVector{1,AbstractSpace}) = spaces[1].size # for `Dimensions.to` has only a single Space
