export AbstractSpace, Space

abstract type AbstractSpace end

# this show function is for printing AbstractDimensions
Base.show(io::IO, svec::SVector{N,<:AbstractSpace}) where {N} = print(io, "[", join(svec, ", "), "]")

struct Space <: AbstractSpace
    size::Int

    function Space(size::Int)
        (size < 1) && throw(DomainError(size, "The size of Space must be positive integer (≥ 1)."))
        return new(size)
    end
end
Base.show(io::IO, s::Space) = print(io, s.size)

# for `prod(::Dimensions)`
Base.:(*)(i::Int, s::Space) = i * s.size
Base.:(*)(s1::Space, s2::Space) = s1.size * s2.size
