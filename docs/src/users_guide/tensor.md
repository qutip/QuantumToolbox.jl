# [Tensor Products and Partial Traces](@id doc:Tensor-products-and-Partial-Traces)

```@setup tensor_products
using QuantumToolbox
```

## [Tensor products](@id doc:Tensor-products)

To describe the states of multipartite quantum systems (such as two coupled qubits, a qubit coupled to an oscillator, etc.) we need to expand the Hilbert space by taking the tensor product of the state vectors for each of the system components. Similarly, the operators acting on the state vectors in the combined Hilbert space (describing the coupled system) are formed by taking the tensor product of the individual operators.

In `QuantumToolbox`, the function [`tensor`](@ref) (or [`kron`](@ref)) is used to accomplish this task. This function takes a collection of [`Ket`](@ref) or [`Operator`](@ref) as argument and returns a composite [`QuantumObject`](@ref) for the combined Hilbert space. The function accepts an arbitrary number of [`QuantumObject`](@ref) as argument. The `type` of returned [`QuantumObject`](@ref) is the same as that of the input(s).

A collection of [`QuantumObject`](@ref):
```@example tensor_products
tensor(sigmax(), sigmax(), sigmax())
```

or a `Vector{QuantumObject}`:

```@example tensor_products
op_list = fill(sigmax(), 3)
tensor(op_list)
```

!!! warning "Beware of type-stability!"
    Please note that `tensor(op_list)` or `kron(op_list)` with `op_list` is a `Vector` is type-instable and can hurt performance. It is recommended to use `tensor(op_list...)` or `kron(op_list...)` instead. See the Section [The Importance of Type-Stability](@ref doc:Type-Stability) for more details.

For example, the state vector describing two qubits in their ground states is formed by taking the tensor product of the two single-qubit ground state vectors:

```@example tensor_products
tensor(basis(2, 0), basis(2, 0))
```

One can generalize to more qubits by adding more component state vectors in the argument list to the [`tensor`](@ref) (or [`kron`](@ref)) function, as illustrated in the following example:

```@example tensor_products
states = QuantumObject[
    normalize(basis(2, 0) + basis(2, 1)),
    normalize(basis(2, 0) + basis(2, 1)),
    basis(2, 0)
]
tensor(states...)
```
This state is slightly more complicated, describing two qubits in a superposition between the up and down states, while the third qubit is in its ground state.

To construct operators that act on an extended Hilbert space of a combined system, we similarly pass a list of operators for each component system to the [`tensor`](@ref) (or [`kron`](@ref)) function. For example, to form the operator that represents the simultaneous action of the ``\hat{\sigma}_x`` operator on two qubits:

```@example tensor_products
tensor(sigmax(), sigmax())
```

To create operators in a combined Hilbert space that only act on a single component, we take the tensor product of the operator acting on the subspace of interest, with the identity operators corresponding to the components that are to be unchanged. For example, the operator that represents ``\hat{\sigma}_z`` on the first qubit in a two-qubit system, while leaving the second qubit unaffected:

```@example tensor_products
tensor(sigmaz(), qeye(2))
```

## Example: Constructing composite Hamiltonians

The [`tensor`](@ref) (or [`kron`](@ref)) function is extensively used when constructing Hamiltonians for composite systems. Here we’ll look at some simple examples.

### Two coupled qubits

First, let’s consider a system of two coupled qubits. Assume that both the qubits have equal energy splitting, and that the qubits are coupled through a ``\hat{\sigma}_x \otimes \hat{\sigma}_x`` interaction with strength ``g = 0.05`` (in units where the bare qubit energy splitting is unity). The Hamiltonian describing this system is:

```@example tensor_products
H = tensor(sigmaz(), qeye(2)) + 
    tensor(qeye(2), sigmaz()) + 
    0.05 * tensor(sigmax(), sigmax())
```

### Three coupled qubits

The two-qubit example is easily generalized to three coupled qubits:

```@example tensor_products
H = tensor(sigmaz(), qeye(2), qeye(2)) + 
    tensor(qeye(2), sigmaz(), qeye(2)) + 
    tensor(qeye(2), qeye(2), sigmaz()) + 
    0.5  * tensor(sigmax(), sigmax(), qeye(2)) + 
    0.25 * tensor(qeye(2), sigmax(), sigmax())
```

### A two-level system coupled to a cavity: The Jaynes-Cummings model

The simplest possible quantum mechanical description for light-matter interaction is encapsulated in the Jaynes-Cummings model, which describes the coupling between a two-level atom and a single-mode electromagnetic field (a cavity mode). Denoting the energy splitting of the atom and cavity ``\omega_a`` and ``\omega_c``, respectively, and the atom-cavity interaction strength ``g``, the Jaynes-Cummings Hamiltonian can be constructed as:

```math
H = \frac{\omega_a}{2}\hat{\sigma}_z + \omega_c \hat{a}^\dagger \hat{a} + g (\hat{a}^\dagger \hat{\sigma}_- + \hat{a} \hat{\sigma}_+)
```

```@example tensor_products
N = 6     # cavity fock space truncation
ωc = 1.25 # frequency of cavity
ωa = 1.0  # frequency of two-level atom
g = 0.75  # interaction strength

a = tensor(qeye(2), destroy(N)) # cavity annihilation operator

# two-level atom operators
σm = tensor(destroy(2), qeye(N))
σz = tensor(sigmaz(), qeye(N))

H = 0.5 * ωa * σz + ωc * a' * a + g * (a' * σm + a * σm')
```

## [Partial trace](@id doc:Partial-trace)

The partial trace is an operation that reduces the dimension of a Hilbert space by eliminating some degrees of freedom by averaging (tracing). In this sense it is therefore the converse of the tensor product. It is useful when one is interested in only a part of a coupled quantum system. For open quantum systems, this typically involves tracing over the environment leaving only the system of interest. In `QuantumToolbox` the function [`ptrace`](@ref) is used to take partial traces. [`ptrace`](@ref) takes one [`QuantumObject`](@ref) as an input, and also one argument `sel`, which marks the component systems that should be kept, and all the other components are traced out. 

Remember that the index of `Julia` starts from `1`, and all the elements in `sel` should be positive `Integer`. Therefore, the type of `sel` can be either `Integer`, `Tuple`, `SVector` ([StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)), or `Vector`.

!!! warning "Beware of type-stability!"
    Although it supports also `Vector` type, it is recommended to use `Tuple` or `SVector` from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) to improve performance. For a brief explanation on the impact of the type of `sel`, see the section [The Importance of Type-Stability](@ref doc:Type-Stability).

For example, the density matrix describing a single qubit obtained from a coupled two-qubit system is obtained via:

```@example tensor_products
ψ = tensor(
    basis(2, 0), 
    basis(2, 1), 
    normalize(basis(2, 0) + basis(2, 1))
)
```

```@example tensor_products
ptrace(ψ, 1) # trace out 2nd and 3rd systems
```

```@example tensor_products
ptrace(ψ, (1, 3)) # trace out 2nd system
```

Note that the partial trace always results in a [`Operator`](@ref) (density matrix), regardless of whether the composite system is a pure state (described by a [`Ket`](@ref)) or a mixed state (described by a [`Operator`](@ref)):

```@example tensor_products
ψ1 = normalize(basis(2, 0) + basis(2, 1))
ψ2 = basis(2, 0)
ψT = tensor(ψ1, ψ2)
```

```@example tensor_products
ptrace(ψT, 1)
```

```@example tensor_products
ρT = tensor(ket2dm(ψ1), ket2dm(ψ1))
```

```@example tensor_products
ptrace(ρT, 1)
```
