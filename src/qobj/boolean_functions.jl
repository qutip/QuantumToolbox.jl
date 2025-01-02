#= 
All boolean functions for checking the data or type in `QuantumObject`
=#

export isket, isbra, isoper, isoperbra, isoperket, issuper
export isunitary

@doc raw"""
    isbra(A)

Checks if the [`QuantumObject`](@ref) `A` is a [`BraQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
isbra(A::QuantumObject) = isbra(typeof(A))
isbra(::Type{QuantumObject{DT,BraQuantumObject,N}}) where {DT,N} = true
isbra(A) = false # default case

@doc raw"""
    isket(A)

Checks if the [`QuantumObject`](@ref) `A` is a [`KetQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
isket(A::QuantumObject) = isket(typeof(A))
isket(::Type{QuantumObject{DT,KetQuantumObject,N}}) where {DT,N} = true
isket(A) = false # default case

@doc raw"""
    isoper(A)

Checks if the [`AbstractQuantumObject`](@ref) `A` is a [`OperatorQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
isoper(A::AbstractQuantumObject) = isoper(typeof(A))
isoper(::Type{<:AbstractQuantumObject{DT,OperatorQuantumObject,N}}) where {DT,N} = true
isoper(A) = false # default case

@doc raw"""
    isoperbra(A)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorBraQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
isoperbra(A::QuantumObject) = isoperbra(typeof(A))
isoperbra(::Type{QuantumObject{DT,OperatorBraQuantumObject,N}}) where {DT,N} = true
isoperbra(A) = false # default case

@doc raw"""
    isoperket(A)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorKetQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
isoperket(A::QuantumObject) = isoperket(typeof(A))
isoperket(::Type{QuantumObject{DT,OperatorKetQuantumObject,N}}) where {DT,N} = true
isoperket(A) = false # default case

@doc raw"""
    issuper(A)

Checks if the [`AbstractQuantumObject`](@ref) `A` is a [`SuperOperatorQuantumObject`](@ref). Default case returns `false` for any other inputs.
"""
issuper(A::AbstractQuantumObject) = issuper(typeof(A))
issuper(::Type{<:AbstractQuantumObject{DT,SuperOperatorQuantumObject,N}}) where {DT,N} = true
issuper(A) = false # default case

@doc raw"""
    ishermitian(A::AbstractQuantumObject)
    isherm(A::AbstractQuantumObject)

Test whether the [`AbstractQuantumObject`](@ref) is Hermitian.

!!! note
    `isherm` is a synonym of `ishermitian`.
"""
LinearAlgebra.ishermitian(A::AbstractQuantumObject) = ishermitian(A.data)

@doc raw"""
    issymmetric(A::AbstractQuantumObject)

Test whether the [`AbstractQuantumObject`](@ref) is symmetric.
"""
LinearAlgebra.issymmetric(A::AbstractQuantumObject) = issymmetric(A.data)

@doc raw"""
    isposdef(A::AbstractQuantumObject)

Test whether the [`AbstractQuantumObject`](@ref) is positive definite (and Hermitian) by trying to perform a Cholesky factorization of `A`.
"""
LinearAlgebra.isposdef(A::AbstractQuantumObject) = isposdef(A.data)

@doc raw"""
    isunitary(U::QuantumObject; kwargs...)

Test whether the [`QuantumObject`](@ref) ``U`` is unitary operator. This function calls `Base.isapprox` to test whether ``U U^\dagger`` is approximately equal to identity operator.

Note that all the keyword arguments will be passed to `Base.isapprox`.
"""
isunitary(U::QuantumObject{<:AbstractArray{T}}; kwargs...) where {T} =
    isoper(U) ? isapprox(U.data * U.data', I(size(U, 1)); kwargs...) : false

@doc raw"""
    SciMLOperators.iscached(A::AbstractQuantumObject)

Test whether the [`AbstractQuantumObject`](@ref) `A` has preallocated caches for inplace evaluations.
"""
SciMLOperators.iscached(A::AbstractQuantumObject) = iscached(A.data)

@doc raw"""
    SciMLOperators.isconstant(A::AbstractQuantumObject)

Test whether the [`AbstractQuantumObject`](@ref) `A` is constant in time. For a [`QuantumObject`](@ref), this function returns `true`, while for a [`QuantumObjectEvolution`](@ref), this function returns `true` if the operator is constant in time.
"""
SciMLOperators.isconstant(A::AbstractQuantumObject) = isconstant(A.data)
