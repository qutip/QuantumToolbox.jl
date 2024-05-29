# [Functions operating on Qobj](@id doc:Functions-operating-on-Qobj)

`QuantumToolbox` also provide functions (methods) that operates on a single [`QuantumObject`](@ref) (taking it as an input):

## Linear algebra

- `A'` or `adjoint(A)`: adjoint (dagger; conjugate transposition)
- `conj`: conjugate
- `transpose`: transpose
- [`partial_transpose`](@ref): partial transpose
- `dot`: dot product
- [`tr`](@ref): trace
- [`ptrace`](@ref): partial trace
- `normalize`: normalization
- `normalize!`: normalization (overwriting input)
- [`expect`](@ref): expectation value
- `sqrt` or [`sqrtm`](@ref): (matrix) square root
- `exp`: (matrix) exponential
- `inv`: (matrix) inverse
- [`svdvals`](@ref): singular values
- [`norm`](@ref): standard vector `p`-norm or [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm
- [`sinm`](@ref): (matrix) sine
- [`cosm`](@ref): (matrix) cosine

## Eigenvalue decomposition

- [`eigenenergies`](@ref): return eigenenergies (values)
- [`eigenstates`](@ref): return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigvals`](@ref): return eigenvalues
- [`eigen`](@ref): using dense eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigsolve`](@ref): using sparse eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigsolve_al`](@ref): using the Arnoldi-Lindblad eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)

## Miscellaneous

- [`tidyup`](@ref): remove small elements from [`QuantumObject`](@ref)
- [`tidyup!`](@ref): remove small elements from [`QuantumObject`](@ref) (overwriting input)
- [`get_data`](@ref): return the data of [`QuantumObject`](@ref)
- [`get_coherence`](@ref): get coherence

## Examples

```@setup Qobj_Function
using QuantumToolbox
```

```@example Qobj_Function
ψ = coherent(5, 1)
```

```@example Qobj_Function
ρ = ψ * ψ'
```

```@example Qobj_Function
tr(ρ)
```

```@example Qobj_Function
norm(ρ)
```

```@example Qobj_Function
sqrtm(ρ)
```

```@example Qobj_Function
normalize(basis(4, 1) + basis(4, 2))
```