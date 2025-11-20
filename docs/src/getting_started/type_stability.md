# [The Importance of Type-Stability](@id doc:Type-Stability)

You are here because you have probably heard about the excellent performance of Julia compared to other common programming languages like Python. One of the reasons is the Just-In-Time (JIT) compiler of Julia, which is able to generate highly optimized machine code. However, the JIT compiler can only do its job if the code type can be inferred. You can also read the [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) section in Julia's documentation for more details. Here, we try to explain it briefly, with a focus on the `QuantumToolbox.jl` package.

!!! note
    This page is not a tutorial on `QuantumToolbox.jl`, but rather a general guide to writing Julia code for simulating quantum systems efficiently. If you don't care about the performance of your code, you can skip this page.

## Basics of type stability

Let's have a look at the following example:

```@setup type-stability
using InteractiveUtils
using QuantumToolbox
```

```@example type-stability
function foo(x)
    if x > 0
        return 1
    else
        return -1.0
    end
end
nothing # hide
```

The function `foo` apparently seems to be innocent. It takes an argument `x` and returns either `1` or `-1.0` depending on the sign of `x`. However, the return type of `foo` is not clear. If `x` is positive, the return type is `Int`, otherwise it is `Float64`. This is a problem for the JIT compiler, because it has to determine the return type of `foo` at runtime. This is called type instability (even though it is a weak form) and may lead to a significant performance penalty. To avoid this, always aim for type-stable code. This means that the return type of a function should be clear from the types of its arguments. We can check the inferred return type of `foo` using the `@code_warntype` macro:

```@example type-stability
@code_warntype foo(1)
```

The key point is to ensure the return type of a function is clear from the types of its arguments. There are several ways to achieve this, and the best approach depends on the specific problem. For example, one can use the same return type:

```@example type-stability
function foo(x)
    if x > 0
        return 1.0
    else
        return -1.0
    end
end
nothing # hide
```

Or you can ensure the return type matches the type of the argument:

```@example type-stability
function foo(x::T) where T
    if x > 0
        return T(1)
    else
        return -T(1)
    end
end
nothing # hide
```

The latter example is very important because it takes advantage of Julia's multiple dispatch, which is one of the most powerful features of the language. Depending on the type `T` of the argument `x`, the Julia compiler generates a specialized version of `foo` that is optimized for this type. If the input type is an `Int64`, the return type is `Int64`, if `x` is a `Float64`, the return type is `Float64`, and so on.

```@example type-stability
@show foo(1)
@show foo(-4.4)
@show foo(1//2)
nothing # hide
```

!!! note
    If you didn't know how to make this function type-stable, it is probably a good idea to read the official Julia documentation, and in particular its [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) section.

## Global variables

Another source of type instability is the use of global variables. In general, it is a good idea to declare global variables as `const` to ensure their type is fixed for the entire program. For example, consider the following function that internally takes a global variable `y`:

```@example type-stability
y = 2.4

function bar(x)
    res = zero(x) # this returns the zero of the same type of x
    for i in 1:1000
        res += y * x
    end
    return res
end
nothing # hide
```

The Julia compiler cannot infer the type of `res` because it depends on the type of `y`, which is a global variable that can change at any time of the program. We can check it using the `@code_warntype` macro:

```@example type-stability
@code_warntype bar(3.2)
```

While in the last example of the `foo` function we got a weak form of type instability, returning a `Union{Int, Float64}`, in this case the return type of `bar` is `Any`, meaning that the compiler doesn't know anything about the return type. Thus, this function has nothing different from a dynamically typed language like Python. We can benchmark the performance of `bar` using the [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) package:

```@example type-stability
using BenchmarkTools

@benchmark bar(3.2)
```

Here we see a lot of memory allocations and low performances in general. To fix this, we can declare a `const` (constant) variable instead:

```@example type-stability
const z = 2.4

function bar(x)
    res = zero(x) # this returns the zero of the same type of x
    for i in 1:1000
        res += z * x
    end
    return res
end

@benchmark bar(3.2)
```

And we can see that the performance has improved significantly. Hence, we highly recommend using global variables as `const`, but only when truly necessary. This choice is problem-dependent, but in the case of `QuantumToolbox.jl`, this can be applied for example in the case of defining the Hilbert space dimensions, static parameters, or the system operators.

Although it is always a good practice to avoid such kind of type instabilities, in the actual implementation of `QuantumToolbox.jl` (where we mainly deal with linear algebra operations), the compiler may perform only a few runtime dispatches, and the performance penalty may be negligible compared to the heavy linear algebra operations.

## Vectors vs Tuples vs StaticArrays

Julia has many ways to represent arrays or lists of general objects. The most common are `Vector`s and `Tuple`s. The former is a dynamic array that can change its size at runtime, while the latter is a fixed-size array that is immutable, and where the type of each element is already known at compile time. For example:

```@example type-stability
v1 = [1, 2, 3] # Vector of Int64
v2 = [1.0 + 2.0im, 3.0 + 4.0im] # Vector of ComplexF64
v3 = [1, "ciao", 3.0] # Vector of Any

t1 = (1, 2, 3) # Tuple of {Int64, Int64, Int64}
t2 = (1.0 + 2.0im, 3.0 + 4.0im) # Tuple of {ComplexF64, ComplexF64}
t3 = (1, "ciao", 3.0) # Tuple of {Int64, String, Float64}

@show typeof(v1)
@show typeof(v2)
@show typeof(v3)
@show typeof(t1)
@show typeof(t2)
@show typeof(t3)
nothing # hide
```

Thus, we highly recommend using `Vector` only when we are sure that it contains elements of the same type, and only when we don't need to know its size at compile time. On the other hand, `Tuple`s are less flexible but more efficient in terms of performance. A third option is to use the `SVector` type from the [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package. This is similar to `Vector`, where the elements should have the same type, but it is fixed-size and immutable. One may ask when it is necessary to know the array size at compile time. A practical example is the case of [`ptrace`](@ref), where it internally reshapes the quantum state into a tensor whose dimensions depend on the number of subsystems. We will see this in more detail in the next section.

## The `QuantumObject` internal structure

Before making a practical example, let's see the internal structure of the [`QuantumObject`](@ref) type. As an example, we consider the case of three qubits, and we study the internal structure of the ``\hat{\sigma}_x^{(2)}`` operator:

```@example type-stability
σx_2 = tensor(qeye(2), sigmax(), qeye(2))
```

and its type is

```@example type-stability
obj_type = typeof(σx_2)
```

This is exactly what the Julia compiler sees: it is a [`QuantumObject`](@ref), composed by a field of type `SparseMatrixCSC{ComplexF64, Int64}` (i.e., the 8x8 matrix containing the Pauli matrix, tensored with the identity matrices of the other two qubits). Then, we can also see that it is a [`Operator`](@ref), with `3` subsystems in total. Hence, just looking at the type of the object, the compiler has all the information it needs to generate a specialized version of the functions.

Let's see more in the details all the internal fields of the [`QuantumObject`](@ref) type:

```@example type-stability
fieldnames(obj_type)
```

```@example type-stability
σx_2.data
```

```@example type-stability
σx_2.type
```


```@example type-stability
σx_2.dims
```

The `dims` field contains the dimensions of the subsystems (in this case, three subsystems with dimension `2` each). We can see that the type of `dims` is `SVector` instead of `Vector`. As we mentioned before, this is very useful in functions like [`ptrace`](@ref). Let's do a simple example of reshaping an operator internally generated from some `dims` input:

```@example type-stability
function reshape_operator_data(dims)
    op = Qobj(randn(hilbert_dimensions_to_size(dims)...), type=Operator(), dims=dims)
    op_dims = op.dims
    op_data = op.data
    return reshape(op_data, vcat(op_dims, op_dims)...)
end

typeof(reshape_operator_data([2, 2, 2]))
```

Which returns a tensor of size `2x2x2x2x2x2`. Let's check the `@code_warntype`:

```@example type-stability
@code_warntype reshape_operator_data([2, 2, 2])
```

We got a `Any` type, because the compiler doesn't know the size of the `dims` vector. We can fix this by using a `Tuple` (or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)):

```@example type-stability
typeof(reshape_operator_data((2, 2, 2)))
```

```@example type-stability
@code_warntype reshape_operator_data((2, 2, 2))
```

Finally, let's look at the benchmarks

```@example type-stability
@benchmark reshape_operator_data($[2, 2, 2])
```

```@example type-stability
@benchmark reshape_operator_data($((2, 2, 2)))
```

Which is an innocuous but huge difference in terms of performance. Hence, we highly recommend using `Tuple` or `SVector` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) when defining the dimensions of a user-defined [`QuantumObject`](@ref).

## The use of `Val` in some `QuantumToolbox.jl` functions

In some functions of `QuantumToolbox.jl`, you may find the use of the [`Val`](https://docs.julialang.org/en/v1/base/base/#Base.Val) type in the arguments. This is a trick to pass a value at compile time, and it is very useful to avoid type instabilities. Let's make a very simple example, where we want to create a Fock state ``|j\rangle`` of a given dimension `N`, and we give the possibility to create it as a sparse or dense vector. At first, we can write the function without using `Val`:

```@example type-stability
using SparseArrays

function my_fock(N::Int, j::Int = 0; sparse::Bool = false)
    if sparse
        array = sparsevec([j + 1], [1.0 + 0im], N)
    else
        array = zeros(ComplexF64, N)
        array[j+1] = 1
    end
    return QuantumObject(array; type = Ket())
end
@show my_fock(2, 1)
@show my_fock(2, 1; sparse = true)
nothing # hide
```

But it is immediately clear that the return type of this function is not clear, because it depends on the value of the `sparse` argument. We can check it using the `@code_warntype` macro:

```@example type-stability
@code_warntype my_fock(2, 1)
```

```@example type-stability
@code_warntype my_fock(2, 1; sparse = true)
```

We can fix this by using the `Val` type, where we enable the multiple dispatch of the function:

```@example type-stability
getVal(::Val{N}) where N = N
function my_fock_good(N::Int, j::Int = 0; sparse::Val = Val(false))
    if getVal(sparse)
        array = zeros(ComplexF64, N)
        array[j+1] = 1
    else
        array = sparsevec([j + 1], [1.0 + 0im], N)
    end
    return QuantumObject(array; type = Ket())
end
@show my_fock_good(2, 1)
@show my_fock_good(2, 1; sparse = Val(true))
nothing # hide
```

And now the return type of the function is clear:

```@example type-stability
@code_warntype my_fock_good(2, 1)
```

```@example type-stability
@code_warntype my_fock_good(2, 1; sparse = Val(true))
```

This is exactly how the current [`fock`](@ref) function is implemented in `QuantumToolbox.jl`. There are many other functions that support this feature, and we highly recommend using it when necessary.

## Conclusions

In this page, we have seen the importance of type stability in Julia, and how to write efficient code in the context of `QuantumToolbox.jl`. We have seen that the internal structure of the [`QuantumObject`](@ref) type is already optimized for the compiler, and we have seen some practical examples of how to write efficient code. We have seen that the use of `Vector` should be avoided when the elements don't have the same type, and that the use of `Tuple` or `SVector` is highly recommended when the size of the array is known at compile time. Finally, we have seen the use of `Val` to pass values at compile time, to avoid type instabilities in some functions.
```

