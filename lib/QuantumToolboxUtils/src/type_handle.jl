makeVal(x::Val{T}) where {T} = x
makeVal(x) = Val(x)

getVal(x::Val{T}) where {T} = T
getVal(x) = x # getVal for any other type

_non_static_array_warning(argname, arg::Tuple{}) =
    throw(ArgumentError("The argument $argname must be a Tuple or a StaticVector of non-zero length."))
_non_static_array_warning(argname, arg::Union{SVector{N, T}, MVector{N, T}, NTuple{N, T}}) where {N, T} = nothing
_non_static_array_warning(argname, arg::AbstractVector{T}) where {T} =
    @warn "The argument $argname should be a Tuple or a StaticVector for better performance. Try to use `$argname = $(Tuple(arg))` instead of `$argname = $arg`. " *
    "Alternatively, you can do `import QuantumToolbox: SVector` " *
    "and use `$argname = SVector(" *
    join(arg, ", ") *
    ")`." maxlog = 1

# lazy tensor warning
for AType in (:AbstractArray, :AbstractSciMLOperator)
    for BType in (:AbstractArray, :AbstractSciMLOperator)
        if AType == BType == :AbstractArray
            @eval begin
                _lazy_tensor_warning(::$AType, ::$BType) = nothing
            end
        else
            @eval begin
                _lazy_tensor_warning(A::$AType, B::$BType) =
                    @warn "using lazy tensor (which can hurt performance) between data types: $(get_typename_wrapper(A)) and $(get_typename_wrapper(B))"
            end
        end
    end
end

get_typename_wrapper(A) = Base.typename(typeof(A)).wrapper

_dense_similar(A::AbstractArray, args...) = similar(A, args...)
_dense_similar(A::AbstractSparseMatrix, args...) = similar(nonzeros(A), args...)

_sparse_similar(A::AbstractArray, args...) = sparse(args...)

# alias of abstract types
const FloatOrComplex = Union{T, Complex{T}} where {T <: AbstractFloat}

# functions for getting Float or Complex element type
_float_type(::AbstractArray{T}) where {T <: Number} = _float_type(T)
_float_type(::AbstractSciMLOperator{T}) where {T <: Number} = _float_type(T)
_float_type(::Type{Int32}) = Float32
_float_type(::Type{Int64}) = Float64
_float_type(::Type{Complex{Int32}}) = Float32
_float_type(::Type{Complex{Int64}}) = Float64
_float_type(::Type{Complex{T}}) where {T <: Real} = T
_float_type(T::Type{<:AbstractFloat}) = T # Allow other untracked Real types, like ForwardDiff.Dual
_complex_float_type(::AbstractArray{T}) where {T <: Number} = _complex_float_type(T)
_complex_float_type(::AbstractSciMLOperator{T}) where {T <: Number} = _complex_float_type(T)
_complex_float_type(::Type{Int32}) = ComplexF32
_complex_float_type(::Type{Int64}) = ComplexF64
_complex_float_type(::Type{Float32}) = ComplexF32
_complex_float_type(::Type{Float64}) = ComplexF64
_complex_float_type(::Type{Complex{Int32}}) = ComplexF32
_complex_float_type(::Type{Complex{Int64}}) = ComplexF64
_complex_float_type(::Type{Complex{Float32}}) = ComplexF32
_complex_float_type(::Type{Complex{Float64}}) = ComplexF64
_complex_float_type(T::Type{<:AbstractFloat}) = Complex{T} # Allow other untracked Complex types, like ForwardDiff.Dual
_complex_float_type(T::Type{<:Complex}) = T       # Allow other untracked Complex types, like ForwardDiff.Dual

_convert_eltype_wordsize(::Type{T}, ::Val{64}) where {T <: Int} = Int64
_convert_eltype_wordsize(::Type{T}, ::Val{32}) where {T <: Int} = Int32
_convert_eltype_wordsize(::Type{T}, ::Val{64}) where {T <: AbstractFloat} = Float64
_convert_eltype_wordsize(::Type{T}, ::Val{32}) where {T <: AbstractFloat} = Float32
_convert_eltype_wordsize(::Type{Complex{T}}, ::Val{64}) where {T <: Union{Int, AbstractFloat}} = ComplexF64
_convert_eltype_wordsize(::Type{Complex{T}}, ::Val{32}) where {T <: Union{Int, AbstractFloat}} = ComplexF32
