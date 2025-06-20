#=
Synonyms of the functions for QuantumObject
=#

export Qobj, QobjEvo, shape, isherm
export trans, dag, matrix_element, unit
export tensor, ⊗
export qeye
export vector_to_operator, operator_to_vector
export sqrtm, logm, expm, sinm, cosm

@doc raw"""
    Qobj(A; kwargs...)

!!! note
    `Qobj` is a synonym for generating [`QuantumObject`](@ref). See the docstring of [`QuantumObject`](@ref) for more details.
"""
const Qobj = QuantumObject # we need the docstring here, otherwise the docstring won't be found because QuantumObject is not a public symbol

@doc raw"""
    QobjEvo(args...; kwargs...)

!!! note
    `QobjEvo` is a synonym for generating [`QuantumObjectEvolution`](@ref). See the docstrings of [`QuantumObjectEvolution`](@ref) for more details.
"""
const QobjEvo = QuantumObjectEvolution # we need the docstring here, otherwise the docstring won't be found because QuantumObjectEvolution is not a public symbol

const shape = size

const isherm = ishermitian

const trans = transpose

const dag = adjoint

@doc raw"""
    matrix_element(i::QuantumObject, A::QuantumObject, j::QuantumObject)

Compute the generalized dot product `dot(i, A*j)` between three [`QuantumObject`](@ref): ``\langle i | \hat{A} | j \rangle``

Note that this function is same as `dot(i, A, j)`

Supports the following inputs:
- `A` is in the type of [`Operator`](@ref), with `i` and `j` are both [`Ket`](@ref).
- `A` is in the type of [`SuperOperator`](@ref), with `i` and `j` are both [`OperatorKet`](@ref)
"""
matrix_element(i, A, j) = dot(i, A, j)

const unit = normalize

const tensor = kron
const ⊗ = kron

const qeye = eye

const vector_to_operator = vec2mat
const operator_to_vector = mat2vec

@doc raw"""
    sqrtm(A::QuantumObject)

Matrix square root of [`Operator`](@ref) type of [`QuantumObject`](@ref)

Note that for other types of [`QuantumObject`](@ref) use `sqrt(A)` instead.
"""
sqrtm(A::QuantumObject{Operator}) = sqrt(A)

@doc raw"""
    logm(A::QuantumObject)

Matrix logarithm of [`QuantumObject`](@ref)

Note that this function is same as `log(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref).
"""
logm(A::QuantumObject{ObjType}) where {ObjType<:Union{Operator,SuperOperator}} = log(A)

@doc raw"""
    expm(A::QuantumObject)

Matrix exponential of [`QuantumObject`](@ref)

Note that this function is same as `exp(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref).
"""
expm(A::QuantumObject{ObjType}) where {ObjType<:Union{Operator,SuperOperator}} = exp(A)

@doc raw"""
    sinm(A::QuantumObject)

Matrix sine of [`QuantumObject`](@ref), defined as

``\sin \left( \hat{A} \right) = \frac{e^{i \hat{A}} - e^{-i \hat{A}}}{2 i}``

Note that this function is same as `sin(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref).
"""
sinm(A::QuantumObject{ObjType}) where {ObjType<:Union{Operator,SuperOperator}} = sin(A)

@doc raw"""
    cosm(A::QuantumObject)

Matrix cosine of [`QuantumObject`](@ref), defined as

``\cos \left( \hat{A} \right) = \frac{e^{i \hat{A}} + e^{-i \hat{A}}}{2}``

Note that this function is same as `cos(A)` and only supports for [`Operator`](@ref) and [`SuperOperator`](@ref).
"""
cosm(A::QuantumObject{ObjType}) where {ObjType<:Union{Operator,SuperOperator}} = cos(A)
