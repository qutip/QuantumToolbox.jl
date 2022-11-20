using LinearAlgebra
using LinearAlgebra: checksquare, BlasFloat, BlasComplex, BlasReal, BlasInt
import LinearAlgebra

abstract type QuantumObjectType end

abstract type BraQuantumObject <: QuantumObjectType end
abstract type KetQuantumObject <: QuantumObjectType end
abstract type OperatorQuantumObject <: QuantumObjectType end
abstract type SuperOperatorQuantumObject <: QuantumObjectType end

mutable struct QuantumObject{MT<:AbstractArray, ObjType<:QuantumObjectType}
    data::MT
    type::Type{ObjType}
    dims::Vector
end

function QuantumObject(A::AbstractVector{T}; type::Type{ObjType} = KetQuantumObject, dims = nothing) where 
        {T, ObjType<:Union{BraQuantumObject, KetQuantumObject}}
        
    dims === nothing ? dims = [length(A)] : nothing
    prod(dims) != length(A) && throw(DimensionMismatch("The dims parameter does not fit the dimension of the Vector."))
    if !(norm(A) ≈ 1)
        @warn "The norm of the input data is not one."
    end
    ObjType <: KetQuantumObject ? QuantumObject(A, type, dims) : QuantumObject(transpose(A), type, dims)
end

function QuantumObject(A::AbstractMatrix{T}; type::Type{ObjType} = OperatorQuantumObject, dims = nothing) where 
        {T, ObjType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}}

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

ket2dm(A::QuantumObject{<:AbstractArray{T}, KetQuantumObject}) where {T} = A * A'

isbra(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = A.type <: BraQuantumObject
isket(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = A.type <: KetQuantumObject
isoper(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = A.type <: OperatorQuantumObject
issuper(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = A.type <: SuperOperatorQuantumObject
Base.size(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = size(A.data)
Base.size(A::QuantumObject{<:AbstractArray{T}, OpType}, inds...) where {T,OpType<:QuantumObjectType} = size(A.data, inds...)
Base.length(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = length(A.data)

SparseArrays.sparse(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T,OpType<:QuantumObjectType} = QuantumObject(sparse(A.data), OpType, A.dims)
SparseArrays.nnz(A::QuantumObject{<:SparseMatrixCSC{T}, OpType}) where {T,OpType<:QuantumObjectType} = nnz(A.data)
SparseArrays.nonzeros(A::QuantumObject{<:SparseMatrixCSC{T}, OpType}) where {T,OpType<:QuantumObjectType} = nonzeros(A.data)
SparseArrays.rowvals(A::QuantumObject{<:SparseMatrixCSC{T}, OpType}) where {T,OpType<:QuantumObjectType} = rowvals(A.data)
SparseArrays.droptol!(A::QuantumObject{<:SparseMatrixCSC{T}, OpType}, tol::Real) where {T,OpType<:QuantumObjectType} = (droptol!(A.data, tol); return A)

Base.isequal(A::QuantumObject{<:AbstractArray{T}, OpType}, B::QuantumObject{<:AbstractArray{T}, OpType}) where 
            {T,OpType<:QuantumObjectType} = isequal(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.isapprox(A::QuantumObject{<:AbstractArray{T}, OpType}, B::QuantumObject{<:AbstractArray{T}, OpType}) where 
            {T,OpType<:QuantumObjectType} = isapprox(A.data, B.data) && isequal(A.type, B.type) && isequal(A.dims, B.dims)
Base.:(==)(A::QuantumObject{<:AbstractArray{T}, OpType}, B::QuantumObject{<:AbstractArray{T}, OpType}) where 
            {T,OpType<:QuantumObjectType} = (A.data == B.data) && (A.type == B.type) && (A.dims == B.dims)
    

LinearAlgebra.Hermitian(A::QuantumObject{<:AbstractArray{T}, OpType}, uplo::Symbol = :U) where 
            {T,OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = 
    QuantumObject(Hermitian(A.data, uplo), A.type, A.dims)

function Base.show(io::IO, ::MIME"text/plain", QO::QuantumObject{<:AbstractArray{T}, OpType}) where 
        {T, OpType<:Union{BraQuantumObject, KetQuantumObject, SuperOperatorQuantumObject}}

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

function Base.show(io::IO, ::MIME"text/plain", QO::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:OperatorQuantumObject}
    op_data = QO.data
    op_dims = QO.dims
    op_type = "Operator"
    println(io, "Quantum Object:   type=", op_type, "   dims=", op_dims, "   size=", size(op_data), "   ishermitian=", ishermitian(op_data))
    show(io, "text/plain", op_data)
end

for op in (:(+), :(-), :(*))
    @eval begin
        function LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T}, OpType}, B::QuantumObject{<:AbstractArray{T}, OpType}) where 
                {T, OpType<:QuantumObjectType}
            
            A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
            QuantumObject($(op)(A.data, B.data), OpType, A.dims)
        end
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = QuantumObject($(op)(A.data), OpType, A.dims)

        LinearAlgebra.$op(n::T1, A::QuantumObject{<:AbstractArray{T2}, OpType}) where {T1<:Number, T2, OpType<:QuantumObjectType} =
            QuantumObject($(op)(n * I, A.data), OpType, A.dims)
        LinearAlgebra.$op(A::QuantumObject{<:AbstractArray{T1}, OpType}, n::T2) where {T1, T2<:Number, OpType<:QuantumObjectType} =
            QuantumObject($(op)(A.data, n * I), OpType, A.dims)
    end
end

function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}, B::QuantumObject{<:AbstractArray{T}, KetQuantumObject}) where {T}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(A.data * B.data, KetQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T}, BraQuantumObject}, B::QuantumObject{<:AbstractArray{T}, OperatorQuantumObject}) where {T}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(A.data * B.data, BraQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T}, KetQuantumObject}, B::QuantumObject{<:AbstractArray{T}, BraQuantumObject}) where {T}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    QuantumObject(sparse(A.data * B.data), OperatorQuantumObject, A.dims)
end
function LinearAlgebra.:(*)(A::QuantumObject{<:AbstractArray{T}, BraQuantumObject}, B::QuantumObject{<:AbstractArray{T}, KetQuantumObject}) where {T}
    A.dims != B.dims && throw(ErrorException("The two operators are not of the same Hilbert dimension."))
    A.data * B.data
end

LinearAlgebra.:(^)(A::QuantumObject{<:AbstractArray{T}, OpType}, n::T1) where {T, T1<:Number, OpType<:QuantumObjectType} = 
    QuantumObject(^(A.data, n), OpType, A.dims)
LinearAlgebra.:(/)(A::QuantumObject{<:AbstractArray{T}, OpType}, n::T1) where {T, T1<:Number, OpType<:QuantumObjectType} = 
    QuantumObject(/(A.data, n), OpType, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}} =
    QuantumObject(adjoint(A.data), OpType, A.dims)
LinearAlgebra.transpose(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}} = 
    QuantumObject(transpose(A.data), OpType, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T}, KetQuantumObject}) where {T} = 
    QuantumObject(adjoint(A.data), BraQuantumObject, A.dims)
LinearAlgebra.adjoint(A::QuantumObject{<:AbstractArray{T}, BraQuantumObject}) where {T} = 
    QuantumObject(adjoint(A.data), KetQuantumObject, A.dims)

LinearAlgebra.tr(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}} = tr(A.data)

LinearAlgebra.norm(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = norm(A.data)
LinearAlgebra.normalize(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = 
    QuantumObject(normalize(A.data), OpType, A.dims)
LinearAlgebra.normalize!(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = 
    (normalize!(A.data); A)
LinearAlgebra.getindex(A::QuantumObject{<:AbstractArray{T}, OpType}, inds...) where {T, OpType<:QuantumObjectType} = getindex(A.data, inds...)
LinearAlgebra.ishermitian(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = ishermitian(A.data)
LinearAlgebra.issymmetric(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = issymmetric(A.data)
LinearAlgebra.isposdef(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = isposdef(A.data)

function LinearAlgebra.kron(A::QuantumObject{<:AbstractArray{T}, OpType}, B::QuantumObject{<:AbstractArray{T}, OpType}) where 
    {T, OpType<:Union{KetQuantumObject, BraQuantumObject, OperatorQuantumObject}}

    QuantumObject(kron(A.data, B.data), OpType, vcat(A.dims, B.dims))
end

LinearAlgebra.exp(A::QuantumObject{<:AbstractArray{T}, OpType}) where {T, OpType<:QuantumObjectType} = 
    QuantumObject(exp(A.data), OpType, A.dims)

LinearAlgebra.triu!(A::QuantumObject{<:AbstractArray{T}, OpType}, k::Int=0) where 
        {T, OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (triu!(A.data, k); A)
LinearAlgebra.tril!(A::QuantumObject{<:AbstractArray{T}, OpType}, k::Int=0) where 
        {T, OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = (tril!(A.data, k); A)
LinearAlgebra.triu(A::QuantumObject{<:AbstractArray{T}, OpType}, k::Int=0) where 
        {T, OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = QuantumObject(triu(A.data, k), OpType, A.dims)
LinearAlgebra.tril(A::QuantumObject{<:AbstractArray{T}, OpType}, k::Int=0) where 
        {T, OpType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject}} = QuantumObject(tril(A.data, k), OpType, A.dims)

function LinearAlgebra.exp(A::SparseMatrixCSC{T,M}; threshold = 1e-14, nonzero_tol = 1e-20) where {T,M}
    rows = checksquare(A) # Throws exception if not square

    mat_norm = norm(A, Inf)
    mat_norm == 0 && return eye(rows).data
    scaling_factor = nextpow(2, mat_norm) # Native routine, faster
    A = A ./ scaling_factor
    delta = 1

    P = spdiagm(0 => ones(T, rows))
    next_term = P
    n = 1

    while delta > threshold
        next_term = T(1 / n) * A * next_term
        droptol!(next_term, nonzero_tol)
        delta = norm(next_term, Inf)
        copyto!(P, P + next_term)
        n = n + 1
    end
    for n = 1:log2(scaling_factor)
        P = P * P;
        if nnz(P) / length(P) < 0.25
            droptol!(P, nonzero_tol)
        end
    end
    P
end

LinearAlgebra.eigen(A::QuantumObject{<:AbstractArray{T}, OpType}) where 
        {T,OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}} = eigen(Array(A.data))
LinearAlgebra.eigvals(A::QuantumObject{<:AbstractArray{T}, OpType}) where 
        {T,OpType<:Union{OperatorQuantumObject, SuperOperatorQuantumObject}} = eigvals(Array(A.data))