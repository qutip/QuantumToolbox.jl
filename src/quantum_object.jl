using LinearAlgebra
using LinearAlgebra: checksquare, BlasFloat, BlasComplex, BlasReal, BlasInt
import LinearAlgebra

abstract type AbstractQuantumObject end
abstract type QuantumObjectType end

@doc raw"""
    BraQuantumObject <: QuantumObjectType

Abstract type representing a bra state ``\bra{\psi}``.
"""
abstract type BraQuantumObject <: QuantumObjectType end

@doc raw"""
    KetQuantumObject <: QuantumObjectType

Abstract type representing a ket state ``\ket{\psi}``.
"""
abstract type KetQuantumObject <: QuantumObjectType end

@doc raw"""
    OperatorQuantumObject <: QuantumObjectType

Abstract type representing an operator ``\hat{O}``.
"""
abstract type OperatorQuantumObject <: QuantumObjectType end

@doc raw"""
    SuperOperatorQuantumObject <: QuantumObjectType

Abstract type representing a super-operator ``\hat{\mathcal{O}}``.
"""
abstract type SuperOperatorQuantumObject <: QuantumObjectType end

@doc raw"""
    mutable struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType}
        data::MT
        type::Type{ObjType}
        dims::Vector{Int}
    end

Julia struct representing any quantum operator.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> a isa QuantumObject
true
```
"""
mutable struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType} <: AbstractQuantumObject
    data::MT
    type::Type{ObjType}
    dims::Vector{Int}
end

function QuantumObject(A::AbstractVector{T}; type::Type{ObjType}=KetQuantumObject, dims=nothing) where
{T,ObjType<:Union{BraQuantumObject,KetQuantumObject}}

    dims === nothing ? dims = [length(A)] : nothing
    prod(dims) != length(A) && throw(DimensionMismatch("The dims parameter does not fit the dimension of the Vector."))
    # if !(norm(A) ≈ 1)
    #     @warn "The norm of the input data is not one."
    # end
    ObjType <: KetQuantumObject ? QuantumObject(A, type, dims) : QuantumObject(transpose(A), type, dims)
end

function QuantumObject(A::AbstractMatrix{T}; type::Type{ObjType}=OperatorQuantumObject, dims=nothing) where
{T,ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}

    n = checksquare(A)
    if type <: SuperOperatorQuantumObject
        sqrt(n) - isqrt(n) != 0 ? throw(DomainError(n, "The dimension of the matrix is not compatible with the SuperOperator type")) : nothing
        dims === nothing ? dims = [isqrt(n)] : nothing
        prod(dims) != isqrt(n) && throw(DimensionMismatch("The dims parameter does not fit the dimension of the Matrix."))
    else
        dims === nothing ? dims = [n] : nothing
        prod(dims) != n && throw(DimensionMismatch("The dims parameter does not fit the dimension of the Matrix."))
    end

    QuantumObject(A, type, dims)
end

function QuantumObject(A::QuantumObject{<:AbstractArray}; type::Type{ObjType}=A.type, dims=A.dims) where
    {ObjType<:QuantumObjectType}

    QuantumObject(A.data, type, dims)
end

@doc raw"""
    Qobj(A::AbstractArray; type::QuantumObjectType, dims::Vector{Int})

Generate `QuantumObject`
"""
Qobj(A; kwargs...) = QuantumObject(A; kwargs...)

@doc raw"""
    ket2dm(ψ::QuantumObject)

Transform the ket state ``\ket{\psi}`` into a pure density matrix ``\hat{\rho} = \dyad{\psi}``.
"""
ket2dm(ψ::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} = ψ * ψ'

ket2dm(ρ::QuantumObject{<:AbstractArray{T},OperatorQuantumObject}) where {T} = ρ

"""
    isbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`BraQuantumObject`](@ref) state.
"""
isbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = A.type <: BraQuantumObject

"""
    isket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`KetQuantumObject`](@ref) state.
"""
isket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = A.type <: KetQuantumObject

"""
    isoper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorQuantumObject`](@ref) state.
"""
isoper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = A.type <: OperatorQuantumObject

"""
    issuper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`SuperOperatorQuantumObject`](@ref) state.
"""
issuper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = A.type <: SuperOperatorQuantumObject

"""
    size(A::QuantumObject)

Returns the size of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.size(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = size(A.data)
Base.size(A::QuantumObject{<:AbstractArray{T},OpType}, inds...) where {T,OpType<:QuantumObjectType} = size(A.data, inds...)

Base.getindex(A::QuantumObject{<:AbstractArray{T},OpType}, inds...) where {T,OpType<:QuantumObjectType} = getindex(A.data, inds...)
Base.setindex!(A::QuantumObject{<:AbstractArray{T},OpType}, val, inds...) where {T,OpType<:QuantumObjectType} = setindex!(A.data, val, inds...)

#    Broadcasting
Base.broadcastable(x::QuantumObject) = x.data
for op in (:(+), :(-), :(*), :(/), :(^))
    @eval begin
        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::QuantumObject)
            return QuantumObject(broadcast($op, x.data, y.data), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::Number)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::Number, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::QuantumObject, y::AbstractArray)
            return QuantumObject(broadcast($op, x.data, y), x.type, x.dims)
        end

        function Base.Broadcast.broadcasted(::typeof($op), x::AbstractArray, y::QuantumObject)
            return QuantumObject(broadcast($op, x, y.data), y.type, y.dims)
        end
    end
end

"""
    length(A::QuantumObject)

Returns the length of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.length(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = length(A.data)

SparseArrays.sparse(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = QuantumObject(sparse(A.data), OpType, A.dims)
SparseArrays.nnz(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject{<:AbstractSparseArray,OpType}, tol::Real) where {OpType<:QuantumObjectType} = (droptol!(A.data, tol); return A)
SparseArrays.dropzeros(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = QuantumObject(dropzeros(A.data), OpType, A.dims)
SparseArrays.dropzeros!(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = (dropzeros!(A.data); return A)

Base.isequal(A::QuantumObject{<:AbstractArray{T},OpType}, B::QuantumObject{<:AbstractArray{T},OpType}) where
{T,OpType<:QuantumObjectType} = isequal(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.isapprox(A::QuantumObject{<:AbstractArray{T},OpType}, B::QuantumObject{<:AbstractArray{T},OpType}) where
{T,OpType<:QuantumObjectType} = isapprox(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.:(==)(A::QuantumObject{<:AbstractArray{T},OpType}, B::QuantumObject{<:AbstractArray{T},OpType}) where
{T,OpType<:QuantumObjectType} = (A.data == B.data) && (A.type == B.type) && (A.dims == B.dims)


LinearAlgebra.Hermitian(A::QuantumObject{<:AbstractArray{T},OpType}, uplo::Symbol=:U) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dims)

function Base.show(io::IO, ::MIME"text/plain", QO::QuantumObject{<:AbstractArray{T},OpType}) where
{T,OpType<:Union{BraQuantumObject,KetQuantumObject,SuperOperatorQuantumObject}}

    op_data = QO.data
    op_dims = QO.dims
    op_type = QO.type
    if op_type <: KetQuantumObject
        op_type = "Ket"
    elseif op_type <: BraQuantumObject
        op_type = "Bra"
    else
        op_type = "SuperOperator"
    end
    println(io, "Quantum Object:   type=", op_type, "   dims=", op_dims, "   size=", size(op_data))
    show(io, "text/plain", op_data)
end

function Base.show(io::IO, ::MIME"text/plain", QO::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:OperatorQuantumObject}
    op_data = QO.data
    op_dims = QO.dims
    op_type = "Operator"
    println(io, "Quantum Object:   type=", op_type, "   dims=", op_dims, "   size=", size(op_data), "   ishermitian=", ishermitian(op_data))
    show(io, "text/plain", op_data)
end

for op in (:(+), :(-), :(*))
    @eval begin
        function LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T1},OpType}, B::QuantumObject{<:AbstractArray{T2},OpType}) where
        {T1,T2,OpType<:QuantumObjectType}

            A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
            QuantumObject($(op)(A.data, B.data), OpType, A.dims)
        end
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = QuantumObject($(op)(A.data), OpType, A.dims)

        LinearAlgebra.$op(n::T1, A::QuantumObject{<:AbstractArray{T2},OpType}) where {T1<:Number,T2,OpType<:QuantumObjectType} =
            QuantumObject($(op)(n * I, A.data), OpType, A.dims)
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T1},OpType}, n::T2) where {T1,T2<:Number,OpType<:QuantumObjectType} =
            QuantumObject($(op)(A.data, n * I), OpType, A.dims)
    end
end

function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject}) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(A.data * B.data, KetQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject}) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(A.data * B.data, BraQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},BraQuantumObject}) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(A.data * B.data, OperatorQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject}) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    A.data * B.data
end

LinearAlgebra.:(^)(A::QuantumObject{<:AbstractArray{T},OpType}, n::T1) where {T,T1<:Number,OpType<:QuantumObjectType} =
    QuantumObject(^(A.data, n), OpType, A.dims)
LinearAlgebra.:(/)(A::QuantumObject{<:AbstractArray{T},OpType}, n::T1) where {T,T1<:Number,OpType<:QuantumObjectType} =
    QuantumObject(/(A.data, n), OpType, A.dims)
LinearAlgebra.dot(A::QuantumObject{<:AbstractArray{T1},OpType},
    B::QuantumObject{<:AbstractArray{T2},OpType}) where {T1<:Number,T2<:Number,OpType<:KetQuantumObject} =
    LinearAlgebra.dot(A.data, B.data)


LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(adjoint(A.data), OpType, A.dims)
LinearAlgebra.transpose(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(transpose(A.data), OpType, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), BraQuantumObject, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},BraQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), KetQuantumObject, A.dims)

LinearAlgebra.inv(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(sparse(inv(Matrix(A.data))), OpType, A.dims)

"""
    tr(A::QuantumObject})

Returns the trace of `A`.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> tr(a' * a)
190.0 + 0.0im
```
"""
LinearAlgebra.tr(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = tr(A.data)

"""
    norm(A::QuantumObject)

Returns the norm of `A`.

# Examples

```
julia> ψ = fock(10, 2)
Quantum Object:   type=Ket   dims=[10]   size=(10,)
10-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> norm(ψ)
1.0
```
"""
LinearAlgebra.norm(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = norm(A.data)
LinearAlgebra.normalize(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(normalize(A.data), OpType, A.dims)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    (normalize!(A.data); A)
LinearAlgebra.ishermitian(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = ishermitian(A.data)
LinearAlgebra.issymmetric(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = issymmetric(A.data)
LinearAlgebra.isposdef(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = isposdef(A.data)

@doc raw"""
    kron(A::QuantumObject, B::QuantumObject)

Returns the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) ``\hat{A} \otimes \hat{B}``.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> kron(a, a)
Quantum Object:   type=Operator   dims=[20, 20]   size=(400, 400)   ishermitian=false
400×400 SparseMatrixCSC{ComplexF64, Int64} with 361 stored entries:
⠀⠀⠘⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠦
```
"""
function LinearAlgebra.kron(A::QuantumObject{<:AbstractArray{T1},OpType}, B::QuantumObject{<:AbstractArray{T2},OpType}) where
{T1,T2,OpType<:Union{KetQuantumObject,BraQuantumObject,OperatorQuantumObject}}

    QuantumObject(kron(A.data, B.data), OpType, vcat(A.dims, B.dims))
end

@doc raw"""
    tensor(A1::QuantumObject, A2::QuantumObject, ...)

Returns the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) ``\hat{A}_1 \otimes \hat{A}_2 \otimes \cdots``.

# Examples

```
julia> x = sigmax()
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 SparseMatrixCSC{ComplexF64, Int64} with 2 stored entries:
     ⋅      1.0+0.0im
 1.0+0.0im      ⋅

julia> x_list = fill(x, 3);

julia> tensor(x_list...)
Quantum Object:   type=Operator   dims=[2, 2, 2]   size=(8, 8)   ishermitian=true
8×8 SparseMatrixCSC{ComplexF64, Int64} with 8 stored entries:
     ⋅          ⋅          ⋅      …      ⋅          ⋅      1.0+0.0im
     ⋅          ⋅          ⋅             ⋅      1.0+0.0im      ⋅    
     ⋅          ⋅          ⋅         1.0+0.0im      ⋅          ⋅    
     ⋅          ⋅          ⋅             ⋅          ⋅          ⋅    
     ⋅          ⋅          ⋅             ⋅          ⋅          ⋅    
     ⋅          ⋅      1.0+0.0im  …      ⋅          ⋅          ⋅    
     ⋅      1.0+0.0im      ⋅             ⋅          ⋅          ⋅    
 1.0+0.0im      ⋅          ⋅             ⋅          ⋅          ⋅
```
"""
tensor(A::QuantumObject...) = kron(A...)

@doc raw"""
    ⊗(A::QuantumObject, B::QuantumObject)

Returns the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) ``\hat{A} \otimes \hat{B}``.

# Examples

```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢

julia> a ⊗ a
Quantum Object:   type=Operator   dims=[20, 20]   size=(400, 400)   ishermitian=false
400×400 SparseMatrixCSC{ComplexF64, Int64} with 361 stored entries:
⠀⠀⠘⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠦
```
"""
⊗(A::QuantumObject, B::QuantumObject) = kron(A, B)

LinearAlgebra.triu!(A::QuantumObject{<:AbstractArray{T},OpType}, k::Integer=0) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (triu!(A.data, k); A)
LinearAlgebra.tril!(A::QuantumObject{<:AbstractArray{T},OpType}, k::Integer=0) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (tril!(A.data, k); A)
LinearAlgebra.triu(A::QuantumObject{<:AbstractArray{T},OpType}, k::Integer=0) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = QuantumObject(triu(A.data, k), OpType, A.dims)
LinearAlgebra.tril(A::QuantumObject{<:AbstractArray{T},OpType}, k::Integer=0) where
{T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = QuantumObject(tril(A.data, k), OpType, A.dims)

LinearAlgebra.lmul!(a::Number, B::QuantumObject{<:AbstractArray}) = (lmul!(a, B.data); B)
LinearAlgebra.rmul!(B::QuantumObject{<:AbstractArray}, a::Number) = (rmul!(B.data, a); B)

@inline LinearAlgebra.mul!(y::AbstractVector{Ty}, A::QuantumObject{<:AbstractMatrix{Ta}}, x, α, β) where {Ty,Ta} = mul!(y, A.data, x, α, β)


LinearAlgebra.sqrt(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(sqrt(A.data), OpType, A.dims)

LinearAlgebra.exp(A::QuantumObject{<:AbstractMatrix{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(dense_to_sparse(exp(A.data)), OpType, A.dims)

LinearAlgebra.exp(A::QuantumObject{<:AbstractSparseMatrix{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(_spexp(A.data), OpType, A.dims)

function _spexp(A::SparseMatrixCSC{T,M}; threshold=1e-14, nonzero_tol=1e-20) where {T,M}
    m = checksquare(A) # Throws exception if not square

    mat_norm = norm(A, Inf)
    mat_norm == 0 && return eye(m).data
    scaling_factor = nextpow(2, mat_norm) # Native routine, faster
    A = A ./ scaling_factor
    delta = 1

    P = spdiagm(0 => ones(T, m))
    next_term = P
    n = 1

    while delta > threshold
        next_term *= A / n
        if nnz(next_term) / length(next_term) > 0.25
            tidyup!(next_term, nonzero_tol)
        end
        delta = norm(next_term, Inf)
        P += next_term
        n += 1
    end
    for n = 1:log2(scaling_factor)
        P = P * P
        if nnz(P) / length(P) > 0.25
            tidyup!(P, nonzero_tol)
        end
    end
    P
end
