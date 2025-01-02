# [Extension for CUDA.jl](@id doc:CUDA)

## Introduction

This is an extension to support `QuantumObject.data` conversion from standard dense and sparse CPU arrays to GPU ([`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl)) arrays.

This extension will be automatically loaded if user imports both `QuantumToolbox.jl` and [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl):

```julia
using QuantumToolbox
using CUDA
using CUDA.CUSPARSE
CUDA.allowscalar(false) # Avoid unexpected scalar indexing
```

We wrapped several functions in `CUDA` and `CUDA.CUSPARSE` in order to not only converting `QuantumObject.data` into GPU arrays, but also changing the element type and word size (`32` and `64`) since some of the GPUs perform better in `32`-bit. The functions are listed as follows (where input `A` is a [`QuantumObject`](@ref)):

- `cu(A; word_size=64)`: return a new [`QuantumObject`](@ref) with `CUDA` arrays and specified `word_size`.
- `CuArray(A)`: If `A.data` is a dense array, return a new [`QuantumObject`](@ref) with `CUDA.CuArray`.
- `CuArray{T}(A)`: If `A.data` is a dense array, return a new [`QuantumObject`](@ref) with `CUDA.CuArray` under element type `T`.
- `CuSparseVector(A)`: If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseVector`.
- `CuSparseVector{T}(A)`: If `A.data` is a sparse vector, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseVector` under element type `T`.
- `CuSparseMatrixCSC(A)`: If `A.data` is a sparse matrix, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseMatrixCSC`.
- `CuSparseMatrixCSC{T}(A)`: If `A.data` is a sparse matrix, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseMatrixCSC` under element type `T`.
- `CuSparseMatrixCSR(A)`: If `A.data` is a sparse matrix, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseMatrixCSR`.
- `CuSparseMatrixCSR{T}(A)`: If `A.data` is a sparse matrix, return a new [`QuantumObject`](@ref) with `CUDA.CUSPARSE.CuSparseMatrixCSR` under element type `T`.

We suggest to convert the arrays from CPU to GPU memory by using the function `cu` because it allows different `data`-types of input [`QuantumObject`](@ref).

Here are some examples:

## Converting dense arrays

```julia
V = fock(2, 0) # CPU dense vector
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 0.0 + 0.0im
```

```julia
cu(V)
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element CuArray{ComplexF64, 1, CUDA.DeviceMemory}:
 1.0 + 0.0im
 0.0 + 0.0im
```

```julia
cu(V; word_size = 32)
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element CuArray{ComplexF32, 1, CUDA.DeviceMemory}:
 1.0 + 0.0im
 0.0 + 0.0im
```

```julia
M = Qobj([1 2; 3 4]) # CPU dense matrix
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false
2×2 Matrix{Int64}:
 1  2
 3  4
```

```julia
cu(M)
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false
2×2 CuArray{Int64, 2, CUDA.DeviceMemory}:
 1  2
 3  4
```

```julia
cu(M; word_size = 32)
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=false
2×2 CuArray{Int32, 2, CUDA.DeviceMemory}:
 1  2
 3  4
```

## Converting sparse arrays

```julia
V = fock(2, 0; sparse=true) # CPU sparse vector
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element SparseVector{ComplexF64, Int64} with 1 stored entry:
  [1]  =  1.0+0.0im
```

```julia
cu(V)
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element CuSparseVector{ComplexF64, Int32} with 1 stored entry:
  [1]  =  1.0+0.0im
```

```julia
cu(V; word_size = 32)
```

```
Quantum Object:   type=Ket   dims=[2]   size=(2,)
2-element CuSparseVector{ComplexF32, Int32} with 1 stored entry:
  [1]  =  1.0+0.0im
```

```julia
M = sigmax() # CPU sparse matrix
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 SparseMatrixCSC{ComplexF64, Int64} with 2 stored entries:
     ⋅      1.0+0.0im
 1.0+0.0im      ⋅    
```

```julia
cu(M)
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 CuSparseMatrixCSC{ComplexF64, Int32} with 2 stored entries:
     ⋅      1.0+0.0im
 1.0+0.0im      ⋅    
```

```julia
cu(M; word_size = 32)
```

```
Quantum Object:   type=Operator   dims=[2]   size=(2, 2)   ishermitian=true
2×2 CuSparseMatrixCSC{ComplexF32, Int32} with 2 stored entries:
     ⋅      1.0+0.0im
 1.0+0.0im      ⋅    
```
