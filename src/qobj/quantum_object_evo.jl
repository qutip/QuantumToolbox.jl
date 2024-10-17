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

where `c_i(p, t)` is a function that depends on the parameters `p` and time `t`, and `\hat{O}_i` are the operators that form the quantum object. The `type` field is the type of the quantum object, and the `dims` field is the dimensions of the quantum object. For more information about `type` and `dims`, see [`QuantumObject`](@ref). For more information about `AbstractSciMLOperator`, see the [SciML](https://docs.sciml.ai/SciMLOperators/stable/) documentation.

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

julia> σm = tensor(qeye(10), sigmam())
Quantum Object:   type=Operator   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 10 stored entries:
⎡⠂⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡀⎦

julia> coef1(p, t) = exp(-1im * t)
coef1 (generic function with 1 method)

julia> coef2(p, t) = sin(t)
coef2 (generic function with 1 method)

julia> op1 = QobjEvo(((a, coef1), (σm, coef2)))
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

julia> op1 = QobjEvo(((a, coef1), (σm, coef2)))
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
"""
struct QuantumObjectEvolution{DT<:AbstractSciMLOperator,ObjType<:QuantumObjectType,N} <:
       AbstractQuantumObject{DT,ObjType,N}
    data::DT
    type::ObjType
    dims::SVector{N,Int}
end

# Make the QuantumObjectEvolution, with the option to pre-multiply by a scalar
function QuantumObjectEvolution(op_func_list::Tuple, α::Union{Nothing,Number} = nothing)
    op, data = _generate_data(op_func_list, α)
    dims = op.dims
    type = op.type

    # Preallocate the SciMLOperator cache using a dense vector as a reference
    v0 = sparse_to_dense(similar(op.data, size(op, 1)))
    data = cache_operator(data, v0)

    return QuantumObjectEvolution(data, type, dims)
end

QuantumObjectEvolution(op::QuantumObject, α::Union{Nothing,Number} = nothing) =
    QuantumObjectEvolution(_make_SciMLOperator(op, α), op.type, op.dims)

function QuantumObjectEvolution(op::QuantumObjectEvolution, α::Union{Nothing,Number} = nothing)
    if α isa Nothing
        return QuantumObjectEvolution(op.data, op.type, op.dims)
    end
    return QuantumObjectEvolution(α * op.data, op.type, op.dims)
end

@generated function _generate_data(op_func_list::Tuple, α)
    op_func_list_types = op_func_list.parameters
    N = length(op_func_list_types)

    dims_expr = ()
    first_op = nothing
    data_expr = :(0)

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

            data_type = op_type.parameters[1]
            dims_expr = (dims_expr..., :(op_func_list[$i][1].dims))
            if i == 1
                first_op = :(op_func_list[$i][1])
            end
        else
            op_type = op_func_type
            (isoper(op_type) || issuper(op_type)) ||
                throw(ArgumentError("The element must be a Operator or SuperOperator."))

            data_type = op_type.parameters[1]
            dims_expr = (dims_expr..., :(op_func_list[$i].dims))

            if i == 1
                first_op = :(op_func_list[$i])
            end
        end
        data_expr = :($data_expr + _make_SciMLOperator(op_func_list[$i], α))
    end

    quote
        dims = tuple($(dims_expr...))

        length(unique(dims)) == 1 || throw(ArgumentError("The dimensions of the operators must be the same."))

        return $first_op, $data_expr
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

function (QO::QuantumObjectEvolution)(p, t)
    # We put 0 in the place of `u` because the time-dependence doesn't depend on the state
    update_coefficients!(QO.data, 0, p, t)
    return QuantumObject(concretize(QO.data), QO.type, QO.dims)
end

(QO::QuantumObjectEvolution)(t) = QO((), t)
