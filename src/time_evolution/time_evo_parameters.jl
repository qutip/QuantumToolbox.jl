# This function should be implemented after Julia v1.12
Base.@constprop :aggressive function _delete_field(a::NamedTuple{an}, field::Symbol) where {an}
    names = Base.diff_names(an, (field,))
    return NamedTuple{names}(a)
end

@doc raw"""
    struct TimeEvolutionParameters

A Julia constructor for handling the parameters of the time evolution of quantum systems.
"""
struct TimeEvolutionParameters{ParT,TE<:AbstractMatrix,PT<:ProgressBar,MCST}
    params::ParT
    expvals::TE
    progr::PT
    mcsolve_params::MCST
end

TimeEvolutionParameters(params, expvals, progr) = TimeEvolutionParameters(params, expvals, progr, nothing)

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
#     return TimeEvolutionParameters(merge(a.params, b), a.expvals, a.progr, a.mcsolve_params)
# end

