# This function should be implemented after Julia v1.12
Base.@constprop :aggressive function _delete_field(a::NamedTuple{an}, field::Symbol) where {an}
    names = Base.diff_names(an, (field,))
    return NamedTuple{names}(a)
end

struct QuantumTimeEvoParameters{TE<:AbstractMatrix,PT<:ProgressBar,ParT}
    expvals::TE
    progr::PT
    params::ParT

    function QuantumTimeEvoParameters(expvals, progr, params)
        _expvals = expvals
        _progr = progr
        _params = params

        # We replace the fields if they are aleady in the `params` struct
        # Then, we remove them from the `params` struct
        if :expvals ∈ fieldnames(typeof(_params))
            _expvals = _params.expvals
            _params = _delete_field(_params, :expvals)
        end
        if :progr ∈ fieldnames(typeof(_params))
            _progr = _params.progr
            _params = _delete_field(_params, :progr)
        end

        return new{typeof(_expvals),typeof(_progr),typeof(_params)}(_expvals, _progr, _params)
    end
end

#=
By defining a custom `getproperty` method for the `QuantumTimeEvoParameters` struct, we can access the fields of `params` directly.
=#
function Base.getproperty(obj::QuantumTimeEvoParameters, field::Symbol)
    if field ∈ fieldnames(typeof(obj))
        getfield(obj, field)
    elseif field ∈ fieldnames(typeof(obj.params))
        getfield(obj.params, field)
    else
        throw(KeyError("Field $field not found in QuantumTimeEvoParameters or params."))
    end
end

#=
It also supports `params` as a `Vector`, so we implement the `getindex` method for the `QuantumTimeEvoParameters` struct.
=#
Base.getindex(obj::QuantumTimeEvoParameters, i::Int) = getindex(obj.params, i)

Base.length(obj::QuantumTimeEvoParameters) = length(obj.params)

function Base.merge(a::QuantumTimeEvoParameters, b::NamedTuple)
    return QuantumTimeEvoParameters(a.expvals, a.progr, merge(a.params, b))
end
