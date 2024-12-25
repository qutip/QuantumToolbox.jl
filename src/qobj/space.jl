export AbstractSpace, Space

abstract type AbstractSpace end

struct Space <: AbstractSpace
    size::Int

    function Space(size::Int)
        (size < 1) && throw(DomainError(size, "The size of Space must be positive integer (≥ 1)."))
        return new(size)
    end
end
Base.string(s::Space) = string(s.size) # this is only used when printing AbstractDimensions

# for `prod(::Dimensions)`
Base.:(*)(i::Int, s::Space) = i * s.size
Base.:(*)(s1::Space, s2::Space) = s1.size * s2.size
