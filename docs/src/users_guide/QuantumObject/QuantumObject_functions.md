# [Functions operating on Qobj](@id doc:Functions-operating-on-Qobj)

`QuantumToolbox` also provide functions (methods) that operates on [`QuantumObject`](@ref).

You can click the function links and see the corresponding docstring for more information.

## Linear algebra and attributes

Here is a table that summarizes all the supported linear algebra functions and attribute functions operating on a given [`QuantumObject`](@ref) `Q`:

| **Description** | **Function call** | **Synonyms** |
|:----------------|:------------------|:-------------|
| conjugate | [`conj(Q)`](@ref conj) | - |
| transpose | [`transpose(Q)`](@ref transpose) | [`trans(Q)`](@ref trans) |
| conjugate transposition | [`adjoint(Q)`](@ref adjoint) | [`Q'`](@ref adjoint), [`dag(Q)`](@ref dag) |
| partial transpose | [`partial_transpose(Q, mask)`](@ref partial_transpose) | - |
| dot product | [`dot(Q1, Q2)`](@ref dot) | - |
| generalized dot product | [`dot(Q1, Q2, Q3)`](@ref dot) | [`matrix_element(Q1, Q2, Q3)`](@ref matrix_element) |
| trace | [`tr(Q)`](@ref tr) | - |
| partial trace | [`ptrace(Q, sel)`](@ref ptrace) | - |
| singular values | [`svdvals(Q)`](@ref svdvals) | - |
| standard vector `p`-norm or [Schatten](https://en.wikipedia.org/wiki/Schatten_norm) `p`-norm | [`norm(Q, p)`](@ref norm) | - |
| normalization | [`normalize(Q, p)`](@ref normalize) | [`unit(Q, p)`](@ref unit) |
| normalization (in-place) | [`normalize!(Q, p)`](@ref normalize!) | - |
| matrix inverse | [`inv(Q)`](@ref inv) | - |
| matrix square root | [`sqrt(Q)`](@ref sqrt) | [`√(Q)`](@ref sqrt), [`sqrtm(Q)`](@ref sqrtm) |
| matrix logarithm | [`log(Q)`](@ref log) | [`logm(Q)`](@ref logm) |
| matrix exponential | [`exp(Q)`](@ref exp) | [`expm(Q)`](@ref expm) |
| matrix sine | [`sin(Q)`](@ref sin) | [`sinm(Q)`](@ref sinm) |
| matrix cosine | [`cos(Q)`](@ref cos) | [`cosm(Q)`](@ref cosm) |
| diagonal elements | [`diag(Q)`](@ref diag) | - |
| projector  | [`proj(Q)`](@ref proj) | - |
| purity | [`purity(Q)`](@ref purity) | - |
| permute | [`permute(Q, order)`](@ref permute) | - |
| remove small elements | [`tidyup(Q, tol)`](@ref tidyup) | - |
| remove small elements (in-place) | [`tidyup!(Q, tol)`](@ref tidyup!) | - |
| get data | [`get_data(Q)`](@ref get_data) | - |
| get coherence | [`get_coherence(Q)`](@ref get_coherence) | - |

## Eigenvalue decomposition

- [`eigenenergies`](@ref): return eigenenergies (eigenvalues)
- [`eigenstates`](@ref): return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigvals`](@ref): return eigenvalues
- [`eigen`](@ref): using dense eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigsolve`](@ref): using sparse eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)
- [`eigsolve_al`](@ref): using the Arnoldi-Lindblad eigen solver and return [`EigsolveResult`](@ref) (contains eigenvalues and eigenvectors)

## Examples

```@setup Qobj_Function
using QuantumToolbox
```

```@example Qobj_Function
ψ = normalize(basis(4, 1) + basis(4, 2))
```

```@example Qobj_Function
ψ'
```

```@example Qobj_Function
ρ = coherent_dm(5, 1)
```

```@example Qobj_Function
diag(ρ)
```

```@example Qobj_Function
get_data(ρ)
```

```@example Qobj_Function
norm(ρ)
```

```@example Qobj_Function
sqrtm(ρ)
```

```@example Qobj_Function
tr(ρ)
```

```@example Qobj_Function
eigenenergies(ρ)
```

```@example Qobj_Function
result = eigenstates(ρ)
```

```@example Qobj_Function
λ, ψ = result
λ # eigenvalues
```

```@example Qobj_Function
ψ # eigenvectors
```

```@example Qobj_Function
λ, ψ, T = result
T # transformation matrix
```