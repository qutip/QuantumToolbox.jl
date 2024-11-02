export QuantumObjectEvolution

@doc raw"""
    struct QuantumObjectEvolution{DT<:AbstractSciMLOperator,ObjType<:QuantumObjectType,N} <: AbstractQuantumObject
        data::DT
        type::ObjType
        dims::SVector{N,Int}
    end

Julia struct representing any time-dependent quantum object. The `data` field is a `AbstractSciMLOperator` object that represents the time-dependent quantum object. It can be seen as

```math
\hat{O}(t) = \sum_{i} c_i(p, t) \hat{O}_i
```

where ``c_i(p, t)`` is a function that depends on the parameters `p` and time `t`, and ``\hat{O}_i`` are the operators that form the quantum object. The `type` field is the type of the quantum object, and the `dims` field is the dimensions of the quantum object. For more information about `type` and `dims`, see [`QuantumObject`](@ref). For more information about `AbstractSciMLOperator`, see the [SciML](https://docs.sciml.ai/SciMLOperators/stable/) documentation.

# Examples
This operator can be initialized in the same way as the QuTiP `QobjEvo` object. For example
```
julia> a = tensor(destroy(10), qeye(2))
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 18 stored entries:
⎡⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⎦

julia> coef1(p, t) = exp(-1im * t)
coef1 (generic function with 1 method)

julia> op = QuantumObjectEvolution(a, coef1)
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=true
ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20)
```

If there are more than 2 operators, we need to put each set of operator and coefficient function into a two-element `Tuple`, and put all these `Tuple`s together in a larger `Tuple`:

```
julia> σm = tensor(qeye(10), sigmam())
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 10 stored entries:
⎡⠂⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡀⎦

julia> coef2(p, t) = sin(t)
coef2 (generic function with 1 method)

julia> op1 = QuantumObjectEvolution(((a, coef1), (σm, coef2)))
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=true
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))
```

We can also concretize the operator at a specific time `t`
```
julia> op1(0.1)
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```

It also supports parameter-dependent time evolution
```
julia> coef1(p, t) = exp(-1im * p.ω1 * t)
coef1 (generic function with 1 method)

julia> coef2(p, t) = sin(p.ω2 * t)
coef2 (generic function with 1 method)

julia> op1 = QuantumObjectEvolution(((a, coef1), (σm, coef2)))
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=true
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))

julia> p = (ω1 = 1.0, ω2 = 0.5)
(ω1 = 1.0, ω2 = 0.5)

julia> op1(p, 0.1)
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```
"""
struct QuantumObjectEvolution{
    DT<:AbstractSciMLOperator,
    ObjType<:Union{OperatorQuantumObject,SuperOperatorQuantumObject},
    N,
} <: AbstractQuantumObject{DT,ObjType,N}
    data::DT
    type::ObjType
    dims::SVector{N,Int}

    function QuantumObjectEvolution(
        data::DT,
        type::ObjType,
        dims,
    ) where {DT<:AbstractSciMLOperator,ObjType<:QuantumObjectType}
        (type == Operator || type == SuperOperator) ||
            throw(ArgumentError("The type $type is not supported for QuantumObjectEvolution."))

        _check_dims(dims)

        _size = _get_size(data)
        _check_QuantumObject(type, dims, _size[1], _size[2])

        N = length(dims)

        return new{DT,ObjType,N}(data, type, SVector{N,Int}(dims))
    end
end

function Base.show(io::IO, QO::QuantumObjectEvolution)
    op_data = QO.data
    println(
        io,
        "Quantum Object Evo.:   type=",
        QO.type,
        "   dims=",
        QO.dims,
        "   size=",
        size(op_data),
        "   ishermitian=",
        ishermitian(op_data),
        "   isconstant=",
        isconstant(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

function QuantumObjectEvolution(data::AbstractSciMLOperator, type::QuantumObjectType, dims::Integer)
    return QuantumObjectEvolution(data, type, SVector{1,Int}(dims))
end

@doc raw"""
    QuantumObjectEvolution(data::AbstractSciMLOperator; type::QuantumObjectType = Operator, dims = nothing)

Generate a [`QuantumObjectEvolution`](@ref) object from a [`SciMLOperator`](https://github.com/SciML/SciMLOperators.jl), in the same way as [`QuantumObject`](@ref) for `AbstractArray` inputs.
"""
function QuantumObjectEvolution(data::AbstractSciMLOperator; type::QuantumObjectType = Operator, dims = nothing)
    _size = _get_size(data)

    if dims isa Nothing
        if type == Operator
            dims = SVector{1,Int}(_size[2])
        elseif type == SuperOperator
            dims = SVector{1,Int}(isqrt(_size[2]))
        end
    end

    return QuantumObjectEvolution(data, type, dims)
end

# Make the QuantumObjectEvolution, with the option to pre-multiply by a scalar
function QuantumObjectEvolution(
    op_func_list::Tuple,
    α::Union{Nothing,Number} = nothing;
    type::Union{Nothing,QuantumObjectType} = nothing,
)
    op, data = _QobjEvo_generate_data(op_func_list, α)
    dims = op.dims
    if type isa Nothing
        type = op.type
    end

    # Preallocate the SciMLOperator cache using a dense vector as a reference
    v0 = sparse_to_dense(similar(op.data, size(op, 1)))
    data = cache_operator(data, v0)

    return QuantumObjectEvolution(data, type, dims)
end

QuantumObjectEvolution(op::QuantumObject, f::Function; type::Union{Nothing,QuantumObjectType} = nothing) =
    QuantumObjectEvolution(((op, f),); type = type)

function QuantumObjectEvolution(
    op::QuantumObject,
    α::Union{Nothing,Number} = nothing;
    type::Union{Nothing,QuantumObjectType} = nothing,
)
    if type isa Nothing
        type = op.type
    end
    return QuantumObjectEvolution(_make_SciMLOperator(op, α), type, op.dims)
end

function QuantumObjectEvolution(
    op::QuantumObjectEvolution,
    α::Union{Nothing,Number} = nothing;
    type::Union{Nothing,QuantumObjectType} = nothing,
)
    if type isa Nothing
        type = op.type
    elseif type != op.type
        throw(
            ArgumentError(
                "The type of the QuantumObjectEvolution object cannot be changed when using another QuantumObjectEvolution object as input.",
            ),
        )
    end
    if α isa Nothing
        return QuantumObjectEvolution(op.data, type, op.dims)
    end
    return QuantumObjectEvolution(α * op.data, type, op.dims)
end

#=
    _QobjEvo_generate_data(op_func_list::Tuple, α; f::Function=identity)

Parse the `op_func_list` and generate the data for the `QuantumObjectEvolution` object. The `op_func_list` is a tuple of tuples or operators. Each element of the tuple can be a tuple with two elements (operator, function) or an operator. The function is used to generate the time-dependent coefficients for the operators. The `α` parameter is used to pre-multiply the operators by a scalar. The `f` parameter is used to pre-applying a function to the operators before converting them to SciML operators. During the parsing, the dimensions of the operators are checked to be the same, and all the constant operators are summed together.

# Arguments
- `op_func_list::Tuple`: A tuple of tuples or operators.
- `α`: A scalar to pre-multiply the operators.
- `f::Function=identity`: A function to pre-apply to the operators.
=#
@generated function _QobjEvo_generate_data(op_func_list::Tuple, α)
    op_func_list_types = op_func_list.parameters
    N = length(op_func_list_types)

    dims_expr = ()
    first_op = nothing
    data_expr = :(0)
    qobj_expr_const = :(0)

    for i in 1:N
        op_func_type = op_func_list_types[i]
        if op_func_type <: Tuple
            length(op_func_type.parameters) == 2 || throw(ArgumentError("The tuple must have two elements."))
            op_type = op_func_type.parameters[1]
            func_type = op_func_type.parameters[2]
            ((isoper(op_type) || issuper(op_type)) && func_type <: Function) || throw(
                ArgumentError(
                    "The first element must be a Operator or SuperOperator, and the second element must be a function.",
                ),
            )

            op = :(op_func_list[$i][1])
            data_type = op_type.parameters[1]
            dims_expr = (dims_expr..., :($op.dims))
            if i == 1
                first_op = :($op)
            end
            data_expr = :($data_expr + _make_SciMLOperator(op_func_list[$i], α))
        else
            op_type = op_func_type
            (isoper(op_type) || issuper(op_type)) ||
                throw(ArgumentError("The element must be a Operator or SuperOperator."))

            data_type = op_type.parameters[1]
            dims_expr = (dims_expr..., :(op_func_list[$i].dims))
            if i == 1
                first_op = :(op_func_list[$i])
            end
            qobj_expr_const = :($qobj_expr_const + op_func_list[$i])
        end
    end

    quote
        dims = tuple($(dims_expr...))

        length(unique(dims)) == 1 || throw(ArgumentError("The dimensions of the operators must be the same."))

        data_expr_const = $qobj_expr_const isa Integer ? $qobj_expr_const : _make_SciMLOperator($qobj_expr_const, α)

        data_expr = data_expr_const + $data_expr

        return $first_op, data_expr
    end
end

function _make_SciMLOperator(op_func::Tuple, α)
    T = eltype(op_func[1])
    update_func = (a, u, p, t) -> op_func[2](p, t)
    if α isa Nothing
        return ScalarOperator(zero(T), update_func) * MatrixOperator(op_func[1].data)
    end
    return ScalarOperator(zero(T), update_func) * MatrixOperator(α * op_func[1].data)
end

function _make_SciMLOperator(op::QuantumObject, α)
    if α isa Nothing
        return MatrixOperator(op.data)
    end
    return MatrixOperator(α * op.data)
end

@doc raw"""
    (A::QuantumObjectEvolution)(ψout, ψin, p, t)

Apply the time-dependent [`QuantumObjectEvolution`](@ref) object `A` to the input state `ψin` at time `t` with parameters `p`. The output state is stored in `ψout`. This function mimics the behavior of a `AbstractSciMLOperator` object.

# Arguments
- `ψout::QuantumObject`: The output state. It must have the same type as `ψin`.
- `ψin::QuantumObject`: The input state. It must be either a [`KetQuantumObject`](@ref) or a [`OperatorKetQuantumObject`](@ref).
- `p`: The parameters of the time-dependent coefficients.
- `t`: The time at which the coefficients are evaluated.

# Returns
- `ψout::QuantumObject`: The output state.

# Examples
```
julia> a = destroy(20)
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 19 stored entries:
⎡⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠈⠢⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⎦

julia> coef1(p, t) = sin(t)
coef1 (generic function with 1 method)

julia> coef2(p, t) = cos(t)
coef2 (generic function with 1 method)

julia> A = QobjEvo(((a, coef1), (a', coef2)))
Quantum Object:   type=Operator   dims=[20]   size=(20, 20)   ishermitian=true
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))

julia> ψ1 = fock(20, 3)
Quantum Object:   type=Ket   dims=[20]   size=(20,)
20-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> ψ2 = zero_ket(20)
Quantum Object:   type=Ket   dims=[20]   size=(20,)
20-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
     ⋮
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im

julia> A(ψ2, ψ1, nothing, 0.1)
20-element Vector{ComplexF64}:
                0.0 + 0.0im
                0.0 + 0.0im
 0.1729165499254989 + 0.0im
                0.0 + 0.0im
 1.9900083305560516 + 0.0im
                    ⋮
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
```
"""
function (A::QuantumObjectEvolution)(
    ψout::QuantumObject{DT1,QobjType},
    ψin::QuantumObject{DT2,QobjType},
    p,
    t,
) where {DT1,DT2,QobjType<:Union{KetQuantumObject,OperatorKetQuantumObject}}
    check_dims(ψout, ψin)
    check_dims(ψout, A)

    if isoper(A) && isoperket(ψin)
        throw(ArgumentError("The input state must be a Ket if the QuantumObjectEvolution object is an Operator."))
    elseif issuper(A) && isket(ψin)
        throw(
            ArgumentError(
                "The input state must be an OperatorKet if the QuantumObjectEvolution object is a SuperOperator.",
            ),
        )
    end

    A.data(ψout.data, ψin.data, p, t)

    return ψout
end

@doc raw"""
    (A::QuantumObjectEvolution)(ψ, p, t)

Apply the time-dependent [`QuantumObjectEvolution`](@ref) object `A` to the input state `ψ` at time `t` with parameters `p`. Out-of-place version of [`(A::QuantumObjectEvolution)(ψout, ψin, p, t)`](@ref). The output state is stored in a new [`QuantumObject`](@ref) object. This function mimics the behavior of a `AbstractSciMLOperator` object.
"""
function (A::QuantumObjectEvolution)(
    ψ::QuantumObject{DT,QobjType},
    p,
    t,
) where {DT,QobjType<:Union{KetQuantumObject,OperatorKetQuantumObject}}
    ψout = QuantumObject(similar(ψ.data), ψ.type, ψ.dims)
    return A(ψout, ψ, p, t)
end

@doc raw"""
    (A::QuantumObjectEvolution)(p, t)

Calculate the time-dependent [`QuantumObjectEvolution`](@ref) object `A` at time `t` with parameters `p`.

# Arguments
- `p`: The parameters of the time-dependent coefficients.
- `t`: The time at which the coefficients are evaluated.

# Returns
- `A::QuantumObject`: The output state.
"""
function (A::QuantumObjectEvolution)(p, t)
    # We put 0 in the place of `u` because the time-dependence doesn't depend on the state
    update_coefficients!(A.data, 0, p, t)
    return QuantumObject(concretize(A.data), A.type, A.dims)
end

(A::QuantumObjectEvolution)(t) = A(nothing, t)

#=
`promote_type` should be applied on types. Here I define `promote_op_type` because it is applied to operators.
=#
promote_op_type(A::QuantumObjectEvolution, B::QuantumObjectEvolution) = get_typename_wrapper(A)
promote_op_type(A::QuantumObjectEvolution, B::QuantumObject) = get_typename_wrapper(A)
promote_op_type(A::QuantumObject, B::QuantumObjectEvolution) = get_typename_wrapper(B)
promote_op_type(A::QuantumObject, B::QuantumObject) = get_typename_wrapper(A)
