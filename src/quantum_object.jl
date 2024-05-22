export AbstractQuantumObject, QuantumObject, Qobj
export QuantumObjectType,
    BraQuantumObject,
    KetQuantumObject,
    OperatorQuantumObject,
    OperatorBraQuantumObject,
    OperatorKetQuantumObject,
    SuperOperatorQuantumObject
export Bra, Ket, Operator, OperatorBra, OperatorKet, SuperOperator

export isket, isbra, isoper, isoperbra, isoperket, issuper, ket2dm
export tensor, ⊗

abstract type AbstractQuantumObject end
abstract type QuantumObjectType end

@doc raw"""
    BraQuantumObject <: QuantumObjectType

Constructor representing a bra state ``\bra{\psi}``.
"""
struct BraQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::BraQuantumObject) = print(io, "Bra")

@doc raw"""
    const Bra = BraQuantumObject()

A constant representing the type of [`BraQuantumObject`](@ref)
"""
const Bra = BraQuantumObject()

@doc raw"""
    KetQuantumObject <: QuantumObjectType

Constructor representing a ket state ``\ket{\psi}``.
"""
struct KetQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::KetQuantumObject) = print(io, "Ket")

@doc raw"""
    const Ket = KetQuantumObject()

A constant representing the type of [`KetQuantumObject`](@ref)
"""
const Ket = KetQuantumObject()

@doc raw"""
    OperatorQuantumObject <: QuantumObjectType

Constructor representing an operator ``\hat{O}``.
"""
struct OperatorQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::OperatorQuantumObject) = print(io, "Operator")

@doc raw"""
    const Operator = OperatorQuantumObject()

A constant representing the type of [`OperatorQuantumObject`](@ref)
"""
const Operator = OperatorQuantumObject()

@doc raw"""
    SuperOperatorQuantumObject <: QuantumObjectType

Constructor representing a super-operator ``\hat{\mathcal{O}}``.
"""
struct SuperOperatorQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::SuperOperatorQuantumObject) = print(io, "SuperOperator")

@doc raw"""
    const SuperOperator = SuperOperatorQuantumObject()

A constant representing the type of [`SuperOperatorQuantumObject`](@ref)
"""
const SuperOperator = SuperOperatorQuantumObject()

@doc raw"""
    OperatorBraQuantumObject <: QuantumObjectType

Constructor representing a bra state in the super-operator formalism ``\langle\langle\rho|``.
"""
struct OperatorBraQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::OperatorBraQuantumObject) = print(io, "OperatorBra")

@doc raw"""
    const OperatorBra = OperatorBraQuantumObject()

A constant representing the type of [`OperatorBraQuantumObject`](@ref)
"""
const OperatorBra = OperatorBraQuantumObject()

@doc raw"""
    OperatorKetQuantumObject <: QuantumObjectType

Constructor representing a ket state in the super-operator formalism ``|\rho\rangle\rangle``.
"""
struct OperatorKetQuantumObject <: QuantumObjectType end
Base.show(io::IO, T::OperatorKetQuantumObject) = print(io, "OperatorKet")

@doc raw"""
    const OperatorKet = OperatorKetQuantumObject()

A constant representing the type of [`OperatorKetQuantumObject`](@ref)
"""
const OperatorKet = OperatorKetQuantumObject()

@doc raw"""
    struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType}
        data::MT
        type::ObjType
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
struct QuantumObject{MT<:AbstractArray,ObjType<:QuantumObjectType} <: AbstractQuantumObject
    data::MT
    type::ObjType
    dims::Vector{Int}
end

function QuantumObject(
    A::AbstractArray{T,N};
    type::ObjType = nothing,
    dims = nothing,
) where {T,N,ObjType<:Union{Nothing,QuantumObjectType}}

    # only accept 1D- and 2D-array
    if N == 1
        Size = (length(A), 1)
    else
        Size = size(A)
        N > 2 ? throw(DomainError(Size, "The dimension of the array is not compatible with Quantum Object")) : nothing
    end

    # decide QuantumObjectType from the size of A
    if type === nothing
        if Size[1] == Size[2]
            type = Operator
        elseif Size[2] == 1
            type = Ket
        elseif Size[1] == 1
            type = Bra
        else
            throw(DomainError(Size, "The dimension of the array is not compatible with Quantum Object"))
        end
    end

    # decide dims from the size of A and the given type
    if dims === nothing
        if (type isa KetQuantumObject) || (type isa OperatorQuantumObject)
            dims = [Size[1]]
        elseif (type isa SuperOperatorQuantumObject) || (type isa OperatorKetQuantumObject)
            dims = [isqrt(Size[1])]
        elseif type isa BraQuantumObject
            dims = [Size[2]]
        elseif type isa OperatorBraQuantumObject
            dims = [isqrt(Size[2])]
        end
    end

    _check_QuantumObject(type, prod(dims), Size[1], Size[2])
    return QuantumObject(A, type, dims)
end

function _check_QuantumObject(type::KetQuantumObject, prod_dims::Int, m::Int, n::Int)
    (n != 1) ? throw(DomainError((m, n), "The dimension of the array is not compatible with Ket type")) : nothing
    return prod_dims != m ? throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) :
           nothing
end

function _check_QuantumObject(type::BraQuantumObject, prod_dims::Int, m::Int, n::Int)
    (m != 1) ? throw(DomainError((m, n), "The dimension of the array is not compatible with Bra type")) : nothing
    return prod_dims != n ? throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) :
           nothing
end

function _check_QuantumObject(type::OperatorQuantumObject, prod_dims::Int, m::Int, n::Int)
    (m != n) ? throw(DomainError((m, n), "The dimension of the array is not compatible with Operator type")) : nothing
    return prod_dims != m ? throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) :
           nothing
end

function _check_QuantumObject(type::SuperOperatorQuantumObject, prod_dims::Int, m::Int, n::Int)
    (m != n) ? throw(DomainError((m, n), "The dimension of the array is not compatible with SuperOperator type")) :
    nothing
    return prod_dims != sqrt(m) ?
           throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) : nothing
end

function _check_QuantumObject(type::OperatorKetQuantumObject, prod_dims::Int, m::Int, n::Int)
    (n != 1) ? throw(DomainError((m, n), "The dimension of the array is not compatible with OperatorKet type")) :
    nothing
    return prod_dims != sqrt(m) ?
           throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) : nothing
end

function _check_QuantumObject(type::OperatorBraQuantumObject, prod_dims::Int, m::Int, n::Int)
    (m != 1) ? throw(DomainError((m, n), "The dimension of the array is not compatible with OperatorBra type")) :
    nothing
    return prod_dims != sqrt(n) ?
           throw(DimensionMismatch("The dims parameter does not fit the dimension of the Array.")) : nothing
end

function QuantumObject(
    A::QuantumObject{<:AbstractArray{T,N}};
    type::ObjType = A.type,
    dims = A.dims,
) where {T,N,ObjType<:QuantumObjectType}
    N == 1 ? Size = (length(A), 1) : Size = size(A)
    _check_QuantumObject(type, prod(dims), Size[1], Size[2])
    return QuantumObject(copy(A.data), type, dims)
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
isbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: BraQuantumObject

"""
    isket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`KetQuantumObject`](@ref) state.
"""
isket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = OpType <: KetQuantumObject

"""
    isoper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorQuantumObject`](@ref) state.
"""
isoper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorQuantumObject

"""
    isoperbra(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorBraQuantumObject`](@ref) state.
"""
isoperbra(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorBraQuantumObject

"""
    isoperket(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`OperatorKetQuantumObject`](@ref) state.
"""
isoperket(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: OperatorKetQuantumObject

"""
    issuper(A::QuantumObject)

Checks if the [`QuantumObject`](@ref) `A` is a [`SuperOperatorQuantumObject`](@ref) state.
"""
issuper(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    OpType <: SuperOperatorQuantumObject

"""
    size(A::QuantumObject)

Returns the size of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.size(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} = size(A.data)
Base.size(A::QuantumObject{<:AbstractArray{T},OpType}, inds...) where {T,OpType<:QuantumObjectType} =
    size(A.data, inds...)

Base.getindex(A::QuantumObject{<:AbstractArray{T},OpType}, inds...) where {T,OpType<:QuantumObjectType} =
    getindex(A.data, inds...)
Base.setindex!(A::QuantumObject{<:AbstractArray{T},OpType}, val, inds...) where {T,OpType<:QuantumObjectType} =
    setindex!(A.data, val, inds...)

"""
    eltype(A::QuantumObject)

Returns the elements type of the matrix or vector corresponding to the [`QuantumObject`](@ref) `A`.
"""
Base.eltype(A::QuantumObject) = eltype(A.data)

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

SparseArrays.sparse(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(sparse(A.data), OpType(), A.dims)
SparseArrays.nnz(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} =
    nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject{<:AbstractSparseArray,OpType}, tol::Real) where {OpType<:QuantumObjectType} =
    (droptol!(A.data, tol); return A)
SparseArrays.dropzeros(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} =
    QuantumObject(dropzeros(A.data), A.type, A.dims)
SparseArrays.dropzeros!(A::QuantumObject{<:AbstractSparseArray,OpType}) where {OpType<:QuantumObjectType} =
    (dropzeros!(A.data); return A)

Base.isequal(
    A::QuantumObject{<:AbstractArray{T},OpType},
    B::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:QuantumObjectType} = isequal(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.isapprox(
    A::QuantumObject{<:AbstractArray{T},OpType},
    B::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:QuantumObjectType} = isapprox(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.:(==)(
    A::QuantumObject{<:AbstractArray{T},OpType},
    B::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:QuantumObjectType} = (A.data == B.data) && (A.type == B.type) && (A.dims == B.dims)

LinearAlgebra.Hermitian(
    A::QuantumObject{<:AbstractArray{T},OpType},
    uplo::Symbol = :U,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dims)

function Base.show(
    io::IO,
    QO::QuantumObject{<:AbstractArray{T},OpType},
) where {
    T,
    OpType<:Union{
        BraQuantumObject,
        KetQuantumObject,
        OperatorBraQuantumObject,
        OperatorKetQuantumObject,
        SuperOperatorQuantumObject,
    },
}
    op_data = QO.data
    println(io, "Quantum Object:   type=", QO.type, "   dims=", QO.dims, "   size=", size(op_data))
    return show(io, MIME("text/plain"), op_data)
end

function Base.show(io::IO, QO::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:OperatorQuantumObject}
    op_data = QO.data
    println(
        io,
        "Quantum Object:   type=",
        QO.type,
        "   dims=",
        QO.dims,
        "   size=",
        size(op_data),
        "   ishermitian=",
        ishermitian(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

for op in (:(+), :(-), :(*))
    @eval begin
        function LinearAlgebra.$op(
            A::QuantumObject{<:AbstractArray{T1},OpType},
            B::QuantumObject{<:AbstractArray{T2},OpType},
        ) where {T1,T2,OpType<:QuantumObjectType}
            A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
            return QuantumObject($(op)(A.data, B.data), OpType(), A.dims)
        end
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
            QuantumObject($(op)(A.data), OpType(), A.dims)

        LinearAlgebra.$op(
            n::T1,
            A::QuantumObject{<:AbstractArray{T2},OpType},
        ) where {T1<:Number,T2,OpType<:QuantumObjectType} = QuantumObject($(op)(n * I, A.data), OpType(), A.dims)
        LinearAlgebra.$op(
            A::QuantumObject{<:AbstractArray{T1},OpType},
            n::T2,
        ) where {T1,T2<:Number,OpType<:QuantumObjectType} = QuantumObject($(op)(A.data, n * I), OpType(), A.dims)
    end
end

function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Ket, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Bra, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},KetQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},BraQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},BraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},KetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(vec2mat(A.data * mat2vec(B.data)), Operator, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorBraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorKetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return A.data * B.data
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},SuperOperatorQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},OperatorKetQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorKet, A.dims)
end
function LinearAlgebra.:(*)(
    A::QuantumObject{<:AbstractArray{T1},OperatorBraQuantumObject},
    B::QuantumObject{<:AbstractArray{T2},SuperOperatorQuantumObject},
) where {T1,T2}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return QuantumObject(A.data * B.data, OperatorBra, A.dims)
end

LinearAlgebra.:(^)(A::QuantumObject{<:AbstractArray{T},OpType}, n::T1) where {T,T1<:Number,OpType<:QuantumObjectType} =
    QuantumObject(^(A.data, n), OpType(), A.dims)
LinearAlgebra.:(/)(A::QuantumObject{<:AbstractArray{T},OpType}, n::T1) where {T,T1<:Number,OpType<:QuantumObjectType} =
    QuantumObject(/(A.data, n), OpType(), A.dims)
function LinearAlgebra.dot(
    A::QuantumObject{<:AbstractArray{T1},OpType},
    B::QuantumObject{<:AbstractArray{T2},OpType},
) where {T1<:Number,T2<:Number,OpType<:Union{KetQuantumObject,OperatorKetQuantumObject}}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    return LinearAlgebra.dot(A.data, B.data)
end

Base.conj(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(conj(A.data), OpType(), A.dims)
LinearAlgebra.adjoint(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(adjoint(A.data), OpType(), A.dims)
LinearAlgebra.transpose(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(transpose(A.data), OpType(), A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},KetQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), Bra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},BraQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), Ket, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},OperatorKetQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), OperatorBra, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T},OperatorBraQuantumObject}) where {T} =
    QuantumObject(adjoint(A.data), OperatorKet, A.dims)

LinearAlgebra.inv(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(sparse(inv(Matrix(A.data))), OpType(), A.dims)

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
LinearAlgebra.tr(
    A::QuantumObject{<:AbstractArray{T},OpType},
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = tr(A.data)

"""
    svdvals(A::QuantumObject)

Return the singular values of a [`QuantumObject`](@ref) in descending order
"""
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractVector}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:DenseMatrix}) = svdvals(A.data)
LinearAlgebra.svdvals(A::QuantumObject{<:AbstractSparseMatrix}) = svdvals(sparse_to_dense(A.data))

"""
    norm(A::QuantumObject, p::Real=2)

If `A` is either [`Ket`](@ref), [`Bra`](@ref), [`OperatorKet`](@ref), or [`OperatorBra`](@ref), returns the standard vector `p`-norm of `A`.
If `A` is either [`Operator`](@ref) or [`SuperOperator`](@ref), returns [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm of `A`.

Note that the default value of `p=2`

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
LinearAlgebra.norm(
    A::QuantumObject{<:AbstractArray{T},OpType},
    p::Real = 2,
) where {T,OpType<:Union{KetQuantumObject,BraQuantumObject,OperatorKetQuantumObject,OperatorBraQuantumObject}} =
    norm(A.data, p)
function LinearAlgebra.norm(
    A::QuantumObject{<:AbstractArray{T},OpType},
    p::Real = 2,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}}
    p == 2.0 && return norm(A.data, 2)
    return norm(svdvals(A), p)
end
LinearAlgebra.normalize(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(normalize(A.data), OpType(), A.dims)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    (normalize!(A.data); A)
LinearAlgebra.ishermitian(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    ishermitian(A.data)
LinearAlgebra.issymmetric(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    issymmetric(A.data)
LinearAlgebra.isposdef(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    isposdef(A.data)

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
function LinearAlgebra.kron(
    A::QuantumObject{<:AbstractArray{T1},OpType},
    B::QuantumObject{<:AbstractArray{T2},OpType},
) where {T1,T2,OpType<:Union{KetQuantumObject,BraQuantumObject,OperatorQuantumObject}}
    return QuantumObject(kron(A.data, B.data), OpType(), vcat(A.dims, B.dims))
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

LinearAlgebra.triu!(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (triu!(A.data, k); A)
LinearAlgebra.tril!(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (tril!(A.data, k); A)
LinearAlgebra.triu(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(triu(A.data, k), OpType(), A.dims)
LinearAlgebra.tril(
    A::QuantumObject{<:AbstractArray{T},OpType},
    k::Integer = 0,
) where {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} =
    QuantumObject(tril(A.data, k), OpType(), A.dims)

LinearAlgebra.lmul!(a::Number, B::QuantumObject{<:AbstractArray}) = (lmul!(a, B.data); B)
LinearAlgebra.rmul!(B::QuantumObject{<:AbstractArray}, a::Number) = (rmul!(B.data, a); B)

@inline LinearAlgebra.mul!(y::AbstractVector{Ty}, A::QuantumObject{<:AbstractMatrix{Ta}}, x, α, β) where {Ty,Ta} =
    mul!(y, A.data, x, α, β)

LinearAlgebra.sqrt(A::QuantumObject{<:AbstractArray{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(sqrt(A.data), OpType(), A.dims)

LinearAlgebra.exp(A::QuantumObject{<:AbstractMatrix{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(dense_to_sparse(exp(A.data)), OpType(), A.dims)

LinearAlgebra.exp(A::QuantumObject{<:AbstractSparseMatrix{T},OpType}) where {T,OpType<:QuantumObjectType} =
    QuantumObject(_spexp(A.data), OpType(), A.dims)

function _spexp(A::SparseMatrixCSC{T,M}; threshold = 1e-14, nonzero_tol = 1e-20) where {T,M}
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
    for n in 1:log2(scaling_factor)
        P = P * P
        if nnz(P) / length(P) > 0.25
            tidyup!(P, nonzero_tol)
        end
    end
    return P
end

# data type conversions
Base.Vector(A::QuantumObject{<:AbstractVector}) = QuantumObject(Vector(A.data), A.type, A.dims)
Base.Vector{T}(A::QuantumObject{<:AbstractVector}) where {T<:Number} = QuantumObject(Vector{T}(A.data), A.type, A.dims)
Base.Matrix(A::QuantumObject{<:AbstractMatrix}) = QuantumObject(Matrix(A.data), A.type, A.dims)
Base.Matrix{T}(A::QuantumObject{<:AbstractMatrix}) where {T<:Number} = QuantumObject(Matrix{T}(A.data), A.type, A.dims)
SparseArrays.SparseVector(A::QuantumObject{<:AbstractVector}) = QuantumObject(SparseVector(A.data), A.type, A.dims)
SparseArrays.SparseVector{T}(A::QuantumObject{<:SparseVector}) where {T<:Number} =
    QuantumObject(SparseVector{T}(A.data), A.type, A.dims)
SparseArrays.SparseMatrixCSC(A::QuantumObject{<:AbstractMatrix}) =
    QuantumObject(SparseMatrixCSC(A.data), A.type, A.dims)
SparseArrays.SparseMatrixCSC{T}(A::QuantumObject{<:SparseMatrixCSC}) where {T<:Number} =
    QuantumObject(SparseMatrixCSC{T}(A.data), A.type, A.dims)
