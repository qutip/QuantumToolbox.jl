# This function should be implemented after Julia v1.12
@doc raw"""
    struct TimeEvolutionParameters

A Julia constructor for handling the parameters of the time evolution of quantum systems.
"""
struct TimeEvolutionParameters{ParT,MCST}
    params::ParT
    mcsolve_params::MCST
end

TimeEvolutionParameters(params) = TimeEvolutionParameters(params, nothing)

#=
By defining a custom `getproperty` method for the `TimeEvolutionParameters` struct, we can access the fields of `params` directly.
=#
function Base.getproperty(obj::TimeEvolutionParameters, field::Symbol)
    if field ∈ fieldnames(typeof(obj))
        getfield(obj, field)
    elseif field ∈ fieldnames(typeof(obj.params))
        getfield(obj.params, field)
    else
        throw(KeyError("Field $field not found in TimeEvolutionParameters or params."))
    end
end

#=
It also supports `params` as a `Vector`, so we implement the `getindex` method for the `TimeEvolutionParameters` struct.
=#
Base.getindex(obj::TimeEvolutionParameters, i::Int) = getindex(obj.params, i)

Base.length(obj::TimeEvolutionParameters) = length(obj.params)

# function Base.merge(a::TimeEvolutionParameters, b::NamedTuple)
#     return TimeEvolutionParameters(merge(a.params, b), a.expvals, a.mcsolve_params)
# end

########## Mark the struct as a SciMLStructure ##########
# The NamedTuple `params` case still doesn't work, and it should be put as a `Vector` instead

isscimlstructure(::TimeEvolutionParameters) = true
# ismutablescimlstructure(::TimeEvolutionParameters{ParT}) where {ParT<:NamedTuple} = false
ismutablescimlstructure(::TimeEvolutionParameters{ParT}) where {ParT<:AbstractVector} = true

hasportion(::Tunable, ::TimeEvolutionParameters) = true

function _vectorize_params(p::TimeEvolutionParameters{ParT}) where {ParT<:NamedTuple}
    buffer = isempty(p.params) ? ComplexF64[] : collect(values(p.params))
    return (buffer, false)
end
_vectorize_params(p::TimeEvolutionParameters{ParT}) where {ParT<:NullParameters} = (ComplexF64[], false)
_vectorize_params(p::TimeEvolutionParameters{ParT}) where {ParT<:AbstractVector} = (p.params, true)

function canonicalize(::Tunable, p::TimeEvolutionParameters)
    buffer, aliases = _vectorize_params(p) # We are assuming that the values have the same type

    # repack takes a new vector of the same length as `buffer`, and constructs
    # a new `TimeEvolutionParameters` object using the values from the new vector for tunables
    # and retaining old values for other parameters. This is exactly what replace does,
    # so we can use that instead.
    repack = let p = p
        repack(newbuffer) = replace(Tunable(), p, newbuffer)
    end
    # the canonicalized vector, the repack function, and a boolean indicating
    # whether the buffer aliases values in the parameter object
    return buffer, repack, aliases
end

function replace(::Tunable, p::TimeEvolutionParameters{ParT}, newbuffer) where {ParT<:NamedTuple}
    @assert length(newbuffer) == length(p.params)
    new_params = NamedTuple{keys(p.params)}(Tuple(newbuffer))
    return TimeEvolutionParameters(new_params, p.mcsolve_params)
end

function replace(::Tunable, p::TimeEvolutionParameters{ParT}, newbuffer) where {ParT<:AbstractVector}
    @assert length(newbuffer) == length(p.params)
    return TimeEvolutionParameters(newbuffer, p.mcsolve_params)
end

function replace!(::Tunable, p::TimeEvolutionParameters{ParT}, newbuffer) where {ParT<:AbstractVector}
    @assert length(newbuffer) == length(p.params)
    copyto!(p.params, newbuffer)
    return p
end
