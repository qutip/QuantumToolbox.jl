export OperatorSum

@doc raw"""
    struct OperatorSum

A constructor to represent a sum of operators ``\sum_i c_i \hat{O}_i`` with a list of coefficients ``c_i`` and a list of operators ``\hat{O}_i``.

This is very useful when we have to update only the coefficients, without allocating memory by performing the sum of the operators.
"""
struct OperatorSum{CT<:AbstractVector{<:Number},OT<:Union{AbstractVector,Tuple}}
    coefficients::CT
    operators::OT
    function OperatorSum(
        coefficients::CT,
        operators::OT,
    ) where {CT<:AbstractVector{<:Number},OT<:Union{AbstractVector,Tuple}}
        length(coefficients) == length(operators) ||
            throw(DimensionMismatch("The number of coefficients must be the same as the number of operators."))
        # Check if all the operators have the same dimensions
        size_1 = size(operators[1])
        mapreduce(x -> size(x) == size_1, &, operators) ||
            throw(DimensionMismatch("All the operators must have the same dimensions."))
        T = promote_type(
            mapreduce(x -> eltype(x), promote_type, operators),
            mapreduce(eltype, promote_type, coefficients),
        )
        coefficients2 = T.(coefficients)
        CT2 = typeof(coefficients2)
        return new{CT2,OT}(coefficients2, operators)
    end
end

Base.size(A::OperatorSum) = size(A.operators[1])
Base.size(A::OperatorSum, i::Int) = size(A.operators[1], i)
Base.size(A::OperatorSum, inds...) = size(A.operators[1], inds...)
Base.length(A::OperatorSum) = length(A.operators[1])
Base.copy(A::OperatorSum) = OperatorSum(copy(A.coefficients), copy(A.operators))
Base.deepcopy(A::OperatorSum) = OperatorSum(deepcopy(A.coefficients), deepcopy(A.operators))

function op_sum_update_coefficients!(A::OperatorSum, coefficients)
    length(A.coefficients) == length(coefficients) ||
        throw(DimensionMismatch("The number of coefficients must be the same as the number of operators."))
    return A.coefficients .= coefficients
end

@inline function LinearAlgebra.mul!(y::AbstractVector{T}, A::OperatorSum, x::AbstractVector, α, β) where {T}
    # Note that β is applied only to the first term
    mul!(y, A.operators[1], x, α * A.coefficients[1], β)
    @inbounds for i in 2:length(A.operators)
        A.coefficients[i] == 0 && continue
        mul!(y, A.operators[i], x, α * A.coefficients[i], 1)
    end
    return y
end

function liouvillian(A::OperatorSum, Id_cache = I(prod(A.operators[1].dims)))
    return OperatorSum(A.coefficients, liouvillian.(A.operators, Ref(Id_cache)))
end
