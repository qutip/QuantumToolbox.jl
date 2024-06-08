#=
Synonyms of the functions for QuantumObject
=#

export Qobj, isherm
export trans, dag, matrix_element, unit
export sqrtm, logm, expm, sinm, cosm

@doc raw"""
    Qobj(A::AbstractArray; type::QuantumObjectType, dims::Vector{Int})

Generate `QuantumObject`

Note that this functions is same as `QuantumObject(A; type=type, dims=dims)`
"""
Qobj(A; kwargs...) = QuantumObject(A; kwargs...)

@doc raw"""
    isherm(A::QuantumObject)

Test whether the [`QuantumObject`](@ref) is Hermitian.

Note that this functions is same as `ishermitian(A)`
"""
isherm(A::QuantumObject{<:AbstractArray{T}}) where {T} = ishermitian(A)

@doc raw"""
    trans(A::QuantumObject)

Lazy matrix transpose of the [`QuantumObject`](@ref).

Note that this function is same as `transpose(A)`
"""
trans(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = transpose(A)

@doc raw"""
    dag(A::QuantumObject)

Lazy adjoint (conjugate transposition) of the [`QuantumObject`](@ref)

Note that this function is same as `adjoint(A)`
"""
dag(A::QuantumObject{<:AbstractArray{T}}) where {T} = adjoint(A)

@doc raw"""
    matrix_element(i::QuantumObject, A::QuantumObject j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between three [`QuantumObject`](@ref): ``\langle i | A | j \rangle``

Note that this function is same as `dot(i, A, j)`

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)
"""
matrix_element(
    i::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},KetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number} = dot(i, A, j)
matrix_element(
    i::QuantumObject{<:AbstractArray{T1},OperatorKetQuantumObject},
    A::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
    j::QuantumObject{<:AbstractArray{T3},OperatorKetQuantumObject},
) where {T1<:Number,T2<:Number,T3<:Number} = dot(i, A, j)

@doc raw"""
    unit(A::QuantumObject, p::Real)

Return normalized [`QuantumObject`](@ref) so that its `p`-norm equals to unity, i.e. `norm(A, p) == 1`.

Support for the following types of [`QuantumObject`](@ref):
- If `A` is [`Ket`](@ref) or [`Bra`](@ref), default `p = 2`
- If `A` is [`Operator`](@ref), default `p = 1`

Note that this function is same as `normalize(A, p)`

Also, see [`norm`](@ref) about its definition for different types of [`QuantumObject`](@ref).
"""
unit(
    A::QuantumObject{<:AbstractArray{T},ObjType},
    p::Real = 2,
) where {T,ObjType<:Union{KetQuantumObject,BraQuantumObject}} = normalize(A, p)
unit(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}, p::Real = 1) where {T} = normalize(A, p)

@doc raw"""
    sqrtm(A::QuantumObject)

Matrix square root of [`Operator`](@ref) type of [`QuantumObject`](@ref)

Note that for other types of [`QuantumObject`](@ref) use `sprt(A)` instead.
"""
sqrtm(A::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = sqrt(A)

@doc raw"""
    logm(A::QuantumObject)

Matrix logarithm of [`QuantumObject`](@ref)

Note that this function is same as `log(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
logm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = log(A)

@doc raw"""
    expm(A::QuantumObject)

Matrix exponential of [`QuantumObject`](@ref)

Note that this function is same as `exp(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
expm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = exp(A)

@doc raw"""
    sinm(A::QuantumObject)

Matrix sine of [`QuantumObject`](@ref), defined as

``\sin \left( \hat{A} \right) = \frac{e^{i \hat{A}} - e^{-i \hat{A}}}{2 i}``

Note that this function is same as `sin(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
sinm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = sin(A)

@doc raw"""
    cosm(A::QuantumObject)

Matrix cosine of [`QuantumObject`](@ref), defined as

``\cos \left( \hat{A} \right) = \frac{e^{i \hat{A}} + e^{-i \hat{A}}}{2}``

Note that this function is same as `cos(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref)
"""
cosm(
    A::QuantumObject{<:AbstractMatrix{T},ObjType},
) where {T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = cos(A)
