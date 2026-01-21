#=
This file defines the QuantumObjectEvolution (QobjEvo) structure.
=#

export QuantumObjectEvolution

@doc raw"""
    struct QuantumObjectEvolution{ObjType<:QuantumObjectType,DimType<:AbstractDimensions,DataType<:AbstractSciMLOperator} <: AbstractQuantumObject{ObjType,DimType,DataType}
        data::DataType
        type::ObjType
        dimensions::DimType
    end

Julia struct representing any time-dependent quantum object. The `data` field is a `AbstractSciMLOperator` object that represents the time-dependent quantum object. It can be seen as

```math
\hat{O}(t) = \sum_{i} c_i(p, t) \hat{O}_i
```

where ``c_i(p, t)`` is a function that depends on the parameters `p` and time `t`, and ``\hat{O}_i`` are the operators that form the quantum object. 

For time-independent cases, see [`QuantumObject`](@ref), and for more information about `AbstractSciMLOperator`, see the [SciML](https://docs.sciml.ai/SciMLOperators/stable/) documentation.

!!! note "`dims` property"
    For a given `H::QuantumObjectEvolution`, `H.dims` or `getproperty(H, :dims)` returns its `dimensions` in the type of integer-vector.

# Examples
This operator can be initialized in the same way as the QuTiP `QobjEvo` object. For example
```jldoctest qobjevo
julia> a = tensor(destroy(10), qeye(2))

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 18 stored entries:
⎡⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⎦

julia> coef1(p, t) = exp(-1im * t)
coef1 (generic function with 1 method)

julia> op = QobjEvo(a, coef1)

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20)
```

If there are more than 2 operators, we need to put each set of operator and coefficient function into a two-element `Tuple`, and put all these `Tuple`s together in a larger `Tuple`:

```jldoctest qobjevo
julia> σm = tensor(qeye(10), sigmam())

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 10 stored entries:
⎡⠂⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡀⎦

julia> coef2(p, t) = sin(t)
coef2 (generic function with 1 method)

julia> op1 = QobjEvo(((a, coef1), (σm, coef2)))

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))
```

We can also concretize the operator at a specific time `t`
```jldoctest qobjevo
julia> op1(0.1)

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```

It also supports parameter-dependent time evolution
```jldoctest qobjevo
julia> coef1(p, t) = exp(-1im * p.ω1 * t)
coef1 (generic function with 1 method)

julia> coef2(p, t) = sin(p.ω2 * t)
coef2 (generic function with 1 method)

julia> op1 = QobjEvo(((a, coef1), (σm, coef2)))

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))

julia> p = (ω1 = 1.0, ω2 = 0.5)
(ω1 = 1.0, ω2 = 0.5)

julia> op1(p, 0.1)

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```
"""
struct QuantumObjectEvolution{
        ObjType <: Union{Operator, SuperOperator},
        DimType <: AbstractDimensions,
        DataType <: AbstractSciMLOperator,
    } <: AbstractQuantumObject{ObjType, DimType, DataType}
    data::DataType
    type::ObjType
    dimensions::DimType

    function QuantumObjectEvolution(data::DT, type, dims) where {DT <: AbstractSciMLOperator}
        ObjType = _check_type(type)
        (type isa Operator || type isa SuperOperator) ||
            throw(ArgumentError("The type $type is not supported for QuantumObjectEvolution."))

        dimensions = _gen_dimensions(type, dims)

        _check_QuantumObject(type, dimensions, size(data, 1), size(data, 2))

        return new{ObjType, typeof(dimensions), DT}(data, type, dimensions)
    end
end

function Base.show(io::IO, QO::QuantumObjectEvolution)
    op_data = QO.data
    println(
        io,
        "\nQuantum Object Evo.:   type=",
        QO.type,
        "   dims=",
        _get_dims_string(QO.dimensions),
        "   size=",
        size(op_data),
        "   ishermitian=",
        ishermitian(op_data),
        "   isconstant=",
        isconstant(op_data),
    )
    return show(io, MIME("text/plain"), op_data)
end

@doc raw"""
    QobjEvo(data::AbstractSciMLOperator; type = Operator(), dims = nothing)
    QuantumObjectEvolution(data::AbstractSciMLOperator; type = Operator(), dims = nothing)

Generate a [`QuantumObjectEvolution`](@ref) object from a [`SciMLOperator`](https://github.com/SciML/SciMLOperators.jl), in the same way as [`QuantumObject`](@ref) for `AbstractArray` inputs.

Note that `QobjEvo` is a synonym of `QuantumObjectEvolution`
"""
function QuantumObjectEvolution(data::AbstractSciMLOperator; type = Operator(), dims = nothing)
    _check_type(type)

    if dims isa Nothing
        if type isa Operator
            dims = ProductDimensions((HilbertSpace(size(data, 1)),), (HilbertSpace(size(data, 2)),))
        elseif type isa SuperOperator
            dims = ProductDimensions((HilbertSpace(isqrt(size(data, 1))),), (HilbertSpace(isqrt(size(data, 2))),))
        end
    end

    return QuantumObjectEvolution(data, type, dims)
end

@doc raw"""
    QobjEvo(op_func_list::Union{Tuple,AbstractQuantumObject}; type=nothing)
    QuantumObjectEvolution(op_func_list::Union{Tuple,AbstractQuantumObject}; type=nothing)

Generate [`QuantumObjectEvolution`](@ref).

# Arguments
- `op_func_list::Union{Tuple,AbstractQuantumObject}`: A tuple of tuples or operators.

!!! warning "Beware of type-stability!"
    Please note that, unlike QuTiP, this function doesn't support `op_func_list` as `Vector` type. This is related to the type-stability issue. See the Section [The Importance of Type-Stability](@ref doc:Type-Stability) for more details.

The `type` parameter is used to specify the type of the [`QuantumObject`](@ref), either `Operator` or `SuperOperator`.

!!! note
    `QobjEvo` is a synonym of `QuantumObjectEvolution`.

# Examples
This operator can be initialized in the same way as the QuTiP `QobjEvo` object. For example
```jldoctest qobjevo
julia> a = tensor(destroy(10), qeye(2))

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 18 stored entries:
⎡⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⎦

julia> σm = tensor(qeye(10), sigmam())

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
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

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))
```

We can also concretize the operator at a specific time `t`
```jldoctest qobjevo
julia> op1(0.1)

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```

It also supports parameter-dependent time evolution
```jldoctest qobjevo
julia> coef1(p, t) = exp(-1im * p.ω1 * t)
coef1 (generic function with 1 method)

julia> coef2(p, t) = sin(p.ω2 * t)
coef2 (generic function with 1 method)

julia> op1 = QobjEvo(((a, coef1), (σm, coef2)))

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))

julia> p = (ω1 = 1.0, ω2 = 0.5)
(ω1 = 1.0, ω2 = 0.5)

julia> op1(p, 0.1)

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 28 stored entries:
⎡⠂⡑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠂⡑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⡑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠂⡑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠂⡑⎦
```
"""
function QuantumObjectEvolution(op_func_list::Tuple; type = nothing)
    op, data = _QobjEvo_generate_data(op_func_list)
    dims = op.dimensions
    _check_type(type)

    if type isa Nothing
        type = op.type
    end

    # Preallocate the SciMLOperator cache using a dense vector as a reference
    v0 = to_dense(similar(op.data, size(op, 1)))
    data = cache_operator(data, v0)

    return QuantumObjectEvolution(data, type, dims)
end

# this is a extra method if user accidentally specify `QuantumObjectEvolution( (op, func) )` or `QuantumObjectEvolution( ((op, func)) )`
QuantumObjectEvolution(op_func::Tuple{<:QuantumObject, <:Function}; type = nothing) =
    QuantumObjectEvolution((op_func,); type = type)

@doc raw"""
    QuantumObjectEvolution(op::QuantumObject, f::Function; type = nothing)
    QobjEvo(op::QuantumObject, f::Function; type = nothing)

Generate [`QuantumObjectEvolution`](@ref).

# Notes
- The `f` parameter is time-dependent coefficient that multiplies the operator. The `type` parameter is used to specify the type of the [`QuantumObject`](@ref), either `Operator` or `SuperOperator`.
- `QobjEvo` is a synonym of `QuantumObjectEvolution`.

# Examples
```jldoctest
julia> a = tensor(destroy(10), qeye(2))

Quantum Object:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=false
20×20 SparseMatrixCSC{ComplexF64, Int64} with 18 stored entries:
⎡⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⎦

julia> coef(p, t) = exp(-1im * t)
coef (generic function with 1 method)

julia> op = QobjEvo(a, coef)

Quantum Object Evo.:   type=Operator()   dims=[10, 2]   size=(20, 20)   ishermitian=true   isconstant=false
ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20)
```
"""
QuantumObjectEvolution(op::QuantumObject, f::Function; type = nothing) = QuantumObjectEvolution(((op, f),); type = type)

function QuantumObjectEvolution(op::QuantumObject; type = nothing)
    _check_type(type)
    if type isa Nothing
        type = op.type
    end
    return QuantumObjectEvolution(_make_SciMLOperator(op), type, op.dimensions)
end

function QuantumObjectEvolution(op::QuantumObjectEvolution; type = nothing)
    _check_type(type)
    if type isa Nothing
        type = op.type
    elseif type != op.type
        throw(
            ArgumentError(
                "The type of the QuantumObjectEvolution object cannot be changed when using another QuantumObjectEvolution object as input.",
            ),
        )
    end
    return QuantumObjectEvolution(op.data, type, op.dimensions)
end

#=
    _QobjEvo_generate_data(op_func_list::Tuple)

Parse the `op_func_list` and generate the data for the `QuantumObjectEvolution` object. The `op_func_list` is a tuple of tuples or operators. Each element of the tuple can be a tuple with two elements (operator, function) or an operator. The function is used to generate the time-dependent coefficients for the operators. During the parsing, the dimensions of the operators are checked to be the same, and all the constant operators are summed together.

# Arguments
- `op_func_list::Tuple`: A tuple of tuples or operators.
=#
function _QobjEvo_generate_data(op_func_list::Tuple)
    first_op = _QobjEvo_get_first_op(op_func_list[1])

    ops_constant = filter(op_func_list) do x
        if x isa QuantumObject
            x.type == first_op.type || throw(ArgumentError("The types of the operators must be the same."))
            x.dimensions == first_op.dimensions ||
                throw(ArgumentError("The dimensions of the operators must be the same."))
            return true
        elseif x isa Tuple
            return false
        else
            throw(ArgumentError("Each element of the tuple must be either a QuantumObject or a Tuple."))
        end
    end

    ops_time_dep = filter(op_func_list) do x
        if x isa Tuple
            _QobjEvo_check_op_func(x)
            op = x[1]
            _QobjEvo_check_op(op)
            op.type == first_op.type || throw(ArgumentError("The types of the operators must be the same."))
            op.dimensions == first_op.dimensions ||
                throw(ArgumentError("The dimensions of the operators must be the same."))
            return true
        elseif x isa QuantumObject
            return false
        else
            throw(ArgumentError("Each element of the tuple must be either a QuantumObject or a Tuple."))
        end
    end

    data_const = isempty(ops_constant) ? zero(eltype(first_op)) : _make_SciMLOperator(sum(ops_constant))
    data_td =
        length(ops_time_dep) == 1 ? _make_SciMLOperator(ops_time_dep[1]) :
        AddedOperator(map(_make_SciMLOperator, ops_time_dep))
    data = data_const + data_td

    return first_op, data
end

function _QobjEvo_check_op(op::AbstractQuantumObject)
    (isoper(op) || issuper(op)) || throw(ArgumentError("The element must be a Operator or SuperOperator."))
    return nothing
end

function _QobjEvo_check_op_func(op_func::Tuple)
    length(op_func) == 2 || throw(ArgumentError("The tuple must have two elements."))
    _QobjEvo_check_op(op_func[1])
    (op_func[2] isa Function) || throw(ArgumentError("The second element must be a function."))
    methods(op_func[2], (Any, Real)) |> length == 0 && throw(
        ArgumentError(
            "The following function must only accept two arguments: `$(nameof(op_func[2]))(p, t)` with t<:Real",
        ),
    )
    return nothing
end

_QobjEvo_get_first_op(op_func_list_1::Union{Tuple, AbstractQuantumObject}) =
if op_func_list_1 isa Tuple
    _QobjEvo_check_op_func(op_func_list_1)
    op = op_func_list_1[1]
    _QobjEvo_check_op(op)
    return op
else
    op = op_func_list_1
    _QobjEvo_check_op(op)
    return op
end

function _make_SciMLOperator(op_func::Tuple)
    op, coef = op_func
    T = eltype(op)
    update_func = (a, u, p, t) -> coef(p, t)
    return ScalarOperator(zero(T), update_func) * _promote_to_scimloperator(op.data)
end
_make_SciMLOperator(op::AbstractQuantumObject) = _promote_to_scimloperator(op.data)

_promote_to_scimloperator(data::AbstractMatrix) = MatrixOperator(data)
_promote_to_scimloperator(data::AbstractSciMLOperator) = data

@doc raw"""
    (A::QuantumObjectEvolution)(ψout, ψin, p, t)

Apply the time-dependent [`QuantumObjectEvolution`](@ref) object `A` to the input state `ψin` at time `t` with parameters `p`. The output state is stored in `ψout`. This function mimics the behavior of a `AbstractSciMLOperator` object.

# Arguments
- `ψout::QuantumObject`: The output state. It must have the same type as `ψin`.
- `ψin::QuantumObject`: The input state. It must be either a [`Ket`](@ref) or a [`OperatorKet`](@ref).
- `p`: The parameters of the time-dependent coefficients.
- `t`: The time at which the coefficients are evaluated.

# Returns
- `ψout::QuantumObject`: The output state.

# Examples
```jldoctest
julia> a = destroy(20)

Quantum Object:   type=Operator()   dims=[20]   size=(20, 20)   ishermitian=false
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

Quantum Object Evo.:   type=Operator()   dims=[20]   size=(20, 20)   ishermitian=true   isconstant=false
(ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20) + ScalarOperator(0.0 + 0.0im) * MatrixOperator(20 × 20))

julia> ψ1 = fock(20, 3);

julia> ψ2 = zero_ket(20);

julia> A(ψ2, ψ1, nothing, 0.1) ≈ A(0.1) * ψ1
true
```
"""
function (A::QuantumObjectEvolution)(
        ψout::QuantumObject{QobjType},
        ψin::QuantumObject{QobjType},
        p,
        t,
    ) where {QobjType <: Union{Ket, OperatorKet}}
    check_mul_dimensions(A, ψin)
    check_dimensions(ψout, ψin)

    if isoper(A) && isoperket(ψin)
        throw(ArgumentError("The input state must be a Ket if the QuantumObjectEvolution object is an Operator."))
    elseif issuper(A) && isket(ψin)
        throw(
            ArgumentError(
                "The input state must be an OperatorKet if the QuantumObjectEvolution object is a SuperOperator.",
            ),
        )
    end

    A.data(ψout.data, ψin.data, nothing, p, t)

    return ψout
end

@doc raw"""
    (A::QuantumObjectEvolution)(ψ, p, t)

Apply the time-dependent [`QuantumObjectEvolution`](@ref) object `A` to the input state `ψ` at time `t` with parameters `p`. Out-of-place version of [`(A::QuantumObjectEvolution)(ψout, ψin, p, t)`](@ref). The output state is stored in a new [`QuantumObject`](@ref) object. This function mimics the behavior of a `AbstractSciMLOperator` object.
"""
function (A::QuantumObjectEvolution)(ψ::QuantumObject{QobjType}, p, t) where {QobjType <: Union{Ket, OperatorKet}}
    ψout = QuantumObject(similar(ψ.data), ψ.type, ψ.dimensions)
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
    return QuantumObject(concretize(A.data), A.type, A.dimensions)
end

(A::QuantumObjectEvolution)(t) = A(nothing, t)

#=
`promote_type` should be applied on types. Here I define `promote_op_type` because it is applied to operators.
=#
promote_op_type(A::QuantumObjectEvolution, B::QuantumObjectEvolution) = get_typename_wrapper(A)
promote_op_type(A::QuantumObjectEvolution, B::QuantumObject) = get_typename_wrapper(A)
promote_op_type(A::QuantumObject, B::QuantumObjectEvolution) = get_typename_wrapper(B)
promote_op_type(A::QuantumObject, B::QuantumObject) = get_typename_wrapper(A)
