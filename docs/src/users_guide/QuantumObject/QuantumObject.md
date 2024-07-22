# [Quantum Objects (Qobj)](@id doc:Qobj)

## Introduction
The key difference between classical and quantum mechanics is the use of operators instead of numbers as variables. Moreover, we need to specify state vectors and their properties. Therefore, in computing the dynamics of quantum systems, we need a data structure that encapsulates the properties of a quantum operator and ket/bra vectors. The quantum object structure, [`QuantumObject`](@ref), accomplishes this using array representation.

[`QuantumObject`](@ref) supports general `Julia` arrays (`Base.AbstractArray`). For example, it can be:
- `Base.Vector` (dense vector)
- `Base.Matrix` (dense matrix)
- `SparseArrays.SparseVector` (sparse vector)
- `SparseArrays.SparseMatrixCSC` (sparse matrix)
- `CUDA.CuArray` (dense GPU vector / matrix)
- `CUDA.CUSPARSE.CuSparseVector` (sparse GPU vector)
- `CUDA.CUSPARSE.CuSparseMatrixCSC` (sparse GPU matrix)
- `CUDA.CUSPARSE.CuSparseMatrixCSR` (sparse GPU matrix)
- and even more ...

We can create a [`QuantumObject`](@ref) with a user defined data set by passing an array of data into the [`QuantumObject`](@ref):

```@setup Qobj
using QuantumToolbox
```

```@example Qobj
QuantumObject([1, 2, 3, 4, 5])
```

We can also use the same function [`Qobj`](@ref) in [`QuTiP` (`Python`)](https://github.com/qutip/qutip):

```@example Qobj
Qobj([1, 2, 3, 4, 5])
```

```@example Qobj
Qobj([1 2 3 4 5])
```

```@example Qobj
Qobj(rand(4, 4))
```

```@example Qobj
Qobj(rand(4, 4), dims = [2, 2])
```

```@example Qobj
Qobj(rand(4, 4), type = SuperOperator)
```

!!! note "Difference between `dims` and `size`"
    Notice that `type`, `dims`, and `size` will change according to the input `data`. Although `dims` and `size` appear to be the same, `dims` keep tracking the dimension of individual Hilbert spaces of a multipartite system, while `size` does not. We refer the reader to the section [tensor products and partial traces](@ref doc:Tensor-products) for more information.

## States and operators

Manually specifying the data for each quantum object is inefficient. Even more so when most objects correspond to commonly used types such as the ladder operators of a harmonic oscillator, the Pauli spin operators for a two-level system, or state vectors such as Fock states. Therefore, `QuantumToolbox` includes predefined functions to construct variety of states and operators (you can click the function links and see the corresponding docstring):

### States
- [`zero_ket`](@ref): zero ket vector
- [`fock`](@ref) or [`basis`](@ref): fock state ket vector
- [`fock_dm`](@ref): density matrix of a fock state
- [`coherent`](@ref): coherent state ket vector 
- [`rand_ket`](@ref): random ket vector
- [`coherent_dm`](@ref): density matrix of a coherent state
- [`thermal_dm`](@ref): density matrix of a thermal state
- [`maximally_mixed_dm`](@ref): density matrix of a maximally mixed state
- [`rand_dm`](@ref): random density matrix
- [`spin_state`](@ref): spin state
- [`spin_coherent`](@ref): coherent spin state
- [`bell_state`](@ref): Bell state
- [`singlet_state`](@ref): two particle singlet state
- [`triplet_states`](@ref): list of the two particle triplet states
- [`w_state`](@ref): `n`-qubit W-state
- [`ghz_state`](@ref): `n`-qudit GHZ state

### Operators
- [`eye`](@ref) or [`qeye`](@ref): identity operator
- [`destroy`](@ref): lowering (destruction) operator
- [`create`](@ref): raising (creation) operator
- [`projection`](@ref): projection operator
- [`displace`](@ref): displacement operator
- [`squeeze`](@ref): single-mode squeeze operator
- [`num`](@ref): bosonic number operator
- [`phase`](@ref): single-mode Pegg-Barnett phase operator
- [`QuantumToolbox.position`](@ref): position operator
- [`QuantumToolbox.momentum`](@ref): momentum operator
- [`sigmax`](@ref): Pauli-``X`` operator
- [`sigmay`](@ref): Pauli-``Y`` operator
- [`sigmaz`](@ref): Pauli-``Z`` operator
- [`sigmap`](@ref): Pauli ladder (``\sigma_+``) operator
- [`sigmam`](@ref): Pauli ladder (``\sigma_-``) operator
- [`jmat`](@ref): high-order Spin-`j` operator
- [`spin_Jx`](@ref): ``S_x`` operator for Spin-`j`
- [`spin_Jy`](@ref): ``S_y`` operator for Spin-`j`
- [`spin_Jz`](@ref): ``S_z`` operator for Spin-`j`
- [`spin_Jm`](@ref): ``S_-`` operator for Spin-`j`
- [`spin_Jp`](@ref): ``S_+`` operator for Spin-`j`
- [`spin_J_set`](@ref): a set of Spin-`j` operators ``(S_x, S_y, S_z)``
- [`fdestroy`](@ref): fermion destruction operator
- [`fcreate`](@ref): fermion creation operator
- [`commutator`](@ref): commutator or anti-commutator
- [`tunneling`](@ref): tunneling operator
- [`qft`](@ref): discrete quantum Fourier transform matrix

As an example, we give the output for a few of these functions:

```@example Qobj
basis(5, 3)
```

```@example Qobj
coherent(5, 0.5 - 0.5im)
```

```@example Qobj
destroy(4)
```

```@example Qobj
sigmaz()
```

## Qobj fields (attributes)

We have seen that a structure [`QuantumObject`](@ref) has several fields (attributes), such as `data`, `type` and `dims`. These can be accessed in the following way:

```@example Qobj
a = destroy(4)
a.data
```

```@example Qobj
a[2, 3] # the indexing in Julia starts from `1`
```

```@example Qobj
a.type
```

```@example Qobj
a.dims
```

In general, the properties of a [`QuantumObject`](@ref) can be retrieved using several functions with inputting [`QuantumObject`](@ref):

```@example Qobj
size(a)
```

```@example Qobj
shape(a) # synonym of size(a)
```

```@example Qobj
length(a)
```

```@example Qobj
eltype(a) # element type
```

```@example Qobj
println(isket(a)) # ket
println(isbra(a)) # bra
println(isoper(a)) # operator
println(isoperket(a)) # operator-ket
println(isoperbra(a)) # operator-bra
println(issuper(a)) # super operator
println(ishermitian(a)) # Hermitian
println(isherm(a)) # synonym of ishermitian(a)
println(issymmetric(a)) # symmetric
println(isposdef(a)) # positive definite (and Hermitian)
```

### `data` conversions

As we mentioned above, `QuantumObject.data` supports general `Julia` arrays. The conversion between different type of `QuantumObject.data` is done using the standard type-conversion for arrays in `Julia`:

```@example Qobj
v_d = basis(2, 0)
```

```@example Qobj
Vector{Int64}(v_d)
```

```@example Qobj
v_s = SparseVector(v_d)
```

```@example Qobj
SparseVector{Float64}(v_s)
```

```@example Qobj
x_s = sigmax()
```

```@example Qobj
SparseMatrixCSC{Int64}(x_s)
```

```@example Qobj
Matrix{Float64}(x_s)
```

To convert between dense and sparse arrays, one can also use [`dense_to_sparse`](@ref) and [`sparse_to_dense`](@ref):

```@example Qobj
x_d = sparse_to_dense(x_s)
```

```@example Qobj
dense_to_sparse(x_d)
```

!!! note "Convert to GPU arrays"
    See [CUDA extension](@ref doc:CUDA) for more details.

`QuantumToolbox` will do conversion when needed to keep everything working in any format. However these conversions could slow down computation and it is recommended to keep to one format family where possible.

## Qobj math

The rules for mathematical operations on [`QuantumObject`](@ref) are similar to the standard scalar, vector, and matrix arithmetic:

```@example Qobj
a = destroy(4)
```

```@example Qobj
a' # synonym of adjoint(a)
```

```@example Qobj
a + 5
```

```@example Qobj
a' * a
```

```@example Qobj
a ^ 3
```

```@example Qobj
x = sigmax()
```

```@example Qobj
x / sqrt(2)
```

```@example Qobj
x ^ 3 == x
```

```@example Qobj
# type \approx + <TAB> to get symbol ≈
x ^ 3 ≈ x 
```

Of course, like matrices, multiplying two objects of incompatible `dims` or `size` throws an error:
```@repl Qobj
a * x
```

Note that there is a special case for multiplying [`Ket`](@ref) and [`Bra`](@ref), which results in outer product ``|u\rangle \otimes \langle v|``:

```@example Qobj
u = Qobj([1, 2, 3])
```

```@example Qobj
v = Qobj([4, 5, 6])
v'
```

```@example Qobj
u * v'
```

Of course, if you switch the order of multiplication, it becomes inner product ``\langle v | u \rangle``:
```@example Qobj
v' * u
```
