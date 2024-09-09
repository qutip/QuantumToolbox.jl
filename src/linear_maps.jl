export AbstractLinearMap

@doc raw"""
    AbstractLinearMap{T, TS}

Represents a general linear map with element type `T` and size `TS`.

## Overview

A **linear map** is a transformation `L` that satisfies:

- **Additivity**: 
    ```math
    L(u + v) = L(u) + L(v)
    ```
- **Homogeneity**: 
    ```math
    L(cu) = cL(u)
    ```

It is typically represented as a matrix with dimensions given by `size`, and this abtract type helps to define this map when the matrix is not explicitly available.

## Methods

- `Base.eltype(A)`: Returns the element type `T`.
- `Base.size(A)`: Returns the size `A.size`.
- `Base.size(A, i)`: Returns the `i`-th dimension.

## Example

As an example, we now define the linear map used in the [`eigsolve_al`](@ref) function for Arnoldi-Lindblad eigenvalue solver:

```julia-repl
struct ArnoldiLindbladIntegratorMap{T,TS,TI} <: AbstractLinearMap{T,TS}
    elty::Type{T}
    size::TS
    integrator::TI
end

function LinearAlgebra.mul!(y::AbstractVector, A::ArnoldiLindbladIntegratorMap, x::AbstractVector)
    reinit!(A.integrator, x)
    solve!(A.integrator)
    return copyto!(y, A.integrator.u)
end
```

where `integrator` is the ODE integrator for the time-evolution. In this way, we can diagonalize this linear map using the [`eigsolve`](@ref) function.
"""
abstract type AbstractLinearMap{T,TS} end

Base.eltype(A::AbstractLinearMap{T}) where {T} = T

Base.size(A::AbstractLinearMap) = A.size
Base.size(A::AbstractLinearMap, i::Int) = A.size[i]
