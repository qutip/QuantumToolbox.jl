export OperatorSum

@doc raw"""
    struct OperatorSum

A structure to represent a sum of operators ``\sum_i c_i \hat{O}_i`` with a list of coefficients ``c_i`` and a list of operators ``\hat{O}_i``.

This is very useful when we have to update only the coefficients, without allocating memory by performing the sum of the operators.
"""
struct OperatorSum{CT<:Vector{<:Number},OT<:AbstractVector} <: AbstractQuantumObject
    coefficients::CT
    operators::OT
    function OperatorSum(coefficients::CT, operators::OT) where {CT<:Vector{<:Number},OT<:AbstractVector}
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
        return new{Vector{T},OT}(coefficients2, operators)
    end
end

Base.size(A::OperatorSum) = size(A.operators[1])
Base.size(A::OperatorSum, inds...) = size(A.operators[1], inds...)
Base.length(A::OperatorSum) = length(A.operators[1])
Base.copy(A::OperatorSum) = OperatorSum(copy(A.coefficients), copy(A.operators))
Base.deepcopy(A::OperatorSum) = OperatorSum(deepcopy(A.coefficients), deepcopy(A.operators))

function update_coefficients!(A::OperatorSum, coefficients)
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
