# [Symmetries and Block Diagonalization](@id doc:Symmetries)

```@setup symmetries
using QuantumToolbox
```

## Introduction

A **symmetry** of a quantum system is an operator ``\hat{\Pi}`` that commutes with its Hamiltonian, ``[\hat{H}, \hat{\Pi}] = 0``. Physically, the quantity represented by ``\hat{\Pi}`` is conserved by the dynamics.

The eigenspaces of ``\hat{\Pi}`` are called **symmetry sectors**. Since ``\hat{H}`` commutes with ``\hat{\Pi}``, it never takes a state out of its sector. Therefore, in a basis made of eigenvectors of ``\hat{\Pi}``, the Hamiltonian is **block diagonal**, with one block for each sector. For two sectors:

```math
\hat{H} = \begin{pmatrix} \hat{H}_{+} & 0 \\ 0 & \hat{H}_{-} \end{pmatrix}.
```

This is useful for two reasons:
- **Cost**: diagonalizing blocks of sizes ``n_1, n_2, \ldots`` costs ``\sum_k n_k^3`` operations instead of ``(\sum_k n_k)^3``.
- **Physics**: every eigenstate belongs to one sector, so we can study each sector separately (e.g., the spectrum of the states with an even number of photons).

This page shows how to find these blocks with [`block_diagonal_form`](@ref), for Hamiltonians and for Liouvillians, and how `QuantumToolbox` describes them with **direct sums** ([`SumSpace`](@ref)).

## Two examples of symmetries

### ``\mathbb{Z}_2``: the parity of two coupled modes

Consider two bosonic modes ``\hat{a}`` and ``\hat{b}`` with Kerr nonlinearities, a two-photon drive on ``\hat{a}``, and a beam splitter between them:

```math
\hat{H} = \Delta (\hat{a}^\dagger \hat{a} + \hat{b}^\dagger \hat{b}) + \frac{U}{2} (\hat{a}^{\dagger 2} \hat{a}^2 + \hat{b}^{\dagger 2} \hat{b}^2) + G (\hat{a}^2 + \hat{a}^{\dagger 2}) + J (\hat{a}^\dagger \hat{b} + \hat{b}^\dagger \hat{a}).
```

The two-photon drive changes ``n_a`` by two, and the beam splitter moves one photon from one mode to the other. No term changes the parity of the total number of photons, so ``\hat{H}`` commutes with the **parity** ``\hat{\Pi} = (-1)^{\hat{a}^\dagger \hat{a} + \hat{b}^\dagger \hat{b}}``. Since ``\hat{\Pi}^2 = 1``, its eigenvalues are ``\pm 1`` and the symmetry group is ``\mathbb{Z}_2 = \{1, \hat{\Pi}\}``. There are two sectors:
- ``V_+``: the Fock states ``|n_a, n_b\rangle`` with ``n_a + n_b`` even;
- ``V_-``: the Fock states ``|n_a, n_b\rangle`` with ``n_a + n_b`` odd.

### ``U(1)``: the number of excitations in the Jaynes-Cummings model

A cavity coupled to a two-level system (a qubit) with the rotating-wave approximation is described by the Jaynes-Cummings Hamiltonian

```math
\hat{H} = \omega_c \hat{a}^\dagger \hat{a} + \omega_q \hat{\sigma}_+ \hat{\sigma}_- + g (\hat{a}^\dagger \hat{\sigma}_- + \hat{a} \hat{\sigma}_+).
```

The interaction exchanges one photon with one excitation of the qubit, so the number of excitations ``\hat{N} = \hat{a}^\dagger \hat{a} + \hat{\sigma}_+ \hat{\sigma}_-`` is conserved. The symmetry operators are ``e^{i \phi \hat{N}}`` for every real ``\phi``, which form the group ``U(1)``. There is one sector ``V_N`` for each number of excitations ``N``: ``V_0`` contains only ``|0, g\rangle``, and ``V_N`` (``N \geq 1``) contains ``|N, g\rangle`` and ``|N - 1, e\rangle``.

## Direct sums

When a space is split into sectors, we write it as a **direct sum**, e.g., ``V = V_+ \oplus V_-``. A vector of ``V_+ \oplus V_-`` is a vector of ``V_+`` followed by a vector of ``V_-``, and an operator on it is a block matrix: the block-diagonal form above is an operator on ``V_+ \oplus V_-`` without off-diagonal blocks.

In `QuantumToolbox`, the direct sum of spaces is a [`SumSpace`](@ref), which can be built with [`direct_sum`](@ref) or its synonym `⊕` (`\oplus` + `<tab>`):

```@example symmetries
S = Space(2) ⊕ Space(3)
```

The direct sum of [`Ket`](@ref)s stacks them, and the direct sum of [`Operator`](@ref)s puts them on the diagonal of a block matrix:

```@example symmetries
ψ = basis(2, 0) ⊕ basis(3, 1)
```

```@example symmetries
A = sigmaz() ⊕ qeye(3)
```

[`direct_component`](@ref) returns a component of a vector, or a block of an operator:

```@example symmetries
direct_component(ψ, 2)
```

The direct sum is associative, so, as for the tensor product, a direct sum of direct sums is flattened:

```@example symmetries
Space(2) ⊕ Space(3) ⊕ Space(4)
```

The `dims` of the quantum objects are printed as `(to, from)`, where:
- a Hilbert space is printed as its size (`3`), and a tensor product with `⊗` (`3⊗3`);
- the Liouville space of the operators that take a ket of ``W`` to a ket of ``V`` (see below) is printed as `ℒ(V ← W)`;
- a direct sum is printed with `⊕`.

The integer form of the `dims` is still given by `A.dims`.

## Block-diagonal form of a Hamiltonian

Let us build the Hamiltonian of the two coupled modes:

```@example symmetries
Na, Nb = 3, 3 # Fock cutoffs
a = destroy(Na) ⊗ qeye(Nb)
b = qeye(Na) ⊗ destroy(Nb)

Δ, U, G, J = -2.0, 1.0, 3.0, 0.5
H = Δ * (a' * a + b' * b) + U / 2 * (a'^2 * a^2 + b'^2 * b^2) + G * (a^2 + a'^2) + J * (a' * b + b' * a)
nothing # hide
```

[`block_diagonal_form`](@ref) finds the blocks from the nonzero elements of the matrix: the basis states connected (directly or through other states) by nonzero elements end up in the same block. We do not need to give it the symmetry operator, but the symmetry must be diagonal in the basis we use, as it is for the parity and the number of excitations in the Fock basis.

```@example symmetries
bdf = block_diagonal_form(H)
bdf.B
```

`bdf.B` is the Hamiltonian in the basis of the blocks: it acts on `5 ⊕ 4`, namely ``V_+`` (5 states, e.g. ``|0, 0\rangle`` and ``|1, 1\rangle``) and ``V_-`` (4 states). `bdf.P` is the permutation from the Fock basis to the basis of the blocks, such that `bdf.B == bdf.P * H * bdf.P'`, and its dims say exactly this: it takes a vector of `3⊗3` to a vector of `5 ⊕ 4`.

```@example symmetries
bdf.P
```

The blocks are in `bdf.blocks`, and their spectra together give the spectrum of `H`:

```@example symmetries
H₊, H₋ = bdf.blocks
sort(vcat(eigenenergies(H₊), eigenenergies(H₋))) ≈ eigenenergies(H)
```

## Open quantum systems

### The Liouvillian

An open quantum system is described by its density matrix ``\hat{\rho}``, which evolves under the Lindblad master equation

```math
\frac{d \hat{\rho}}{dt} = \mathcal{L} \hat{\rho} = -i [\hat{H}, \hat{\rho}] + \sum_k \left( \hat{C}_k \hat{\rho} \hat{C}_k^\dagger - \frac{1}{2} \{\hat{C}_k^\dagger \hat{C}_k, \hat{\rho}\} \right),
```

where ``\hat{C}_k`` are the jump (or collapse) operators, e.g. ``\hat{C} = \sqrt{\kappa} \hat{a}`` for the loss of photons. The **Liouvillian** ``\mathcal{L}`` is a linear map acting on operators, a [`SuperOperator`](@ref) built by [`liouvillian`](@ref). To write it as a matrix, the density matrix is turned into a vector by stacking its columns ([`mat2vec`](@ref)): this vector belongs to the **Liouville space** ``\mathcal{L}(V \leftarrow V)`` of the operators on ``V`` (a [`LiouvilleSpace`](@ref)), whose dimension is ``\dim(V)^2``.

### Strong and weak symmetries

With the sectors ``V_+`` and ``V_-``, the density matrix is made of four blocks:

```math
\hat{\rho} = \begin{pmatrix} \hat{\rho}_{++} & \hat{\rho}_{+-} \\ \hat{\rho}_{-+} & \hat{\rho}_{--} \end{pmatrix},
```

where ``\hat{\rho}_{+-}`` takes a state of ``V_-`` to a state of ``V_+``. The Hamiltonian part ``-i[\hat{H}, \hat{\rho}]`` maps each block to itself, while the jumps depend on how the jump operators behave under the symmetry:

- **Strong symmetry**: every jump operator commutes with ``\hat{\Pi}``, e.g. the loss of photons in pairs, ``\hat{C} = \hat{a}^2``. Then each block of ``\hat{\rho}`` evolves alone, and the Liouvillian has **four** blocks. Both ``\hat{\rho}_{++}`` and ``\hat{\rho}_{--}`` keep their trace, so there is a steady state in each of them: the final state depends on the initial parity.
- **Weak symmetry**: every jump operator only picks a phase, ``\hat{\Pi} \hat{C}_k \hat{\Pi}^\dagger = e^{i \theta_k} \hat{C}_k``, e.g. the loss of single photons, ``\hat{\Pi} \hat{a} \hat{\Pi}^\dagger = -\hat{a}``. The jump ``\hat{a} \hat{\rho} \hat{a}^\dagger`` moves ``\hat{\rho}_{--}`` to ``\hat{\rho}_{++}`` and vice versa, and ``\hat{\rho}_{+-}`` to ``\hat{\rho}_{-+}``. Then the Liouvillian has **two** blocks: the even one, made of ``\hat{\rho}_{++}`` and ``\hat{\rho}_{--}``, which contains the steady state, and the odd one, made of the coherences ``\hat{\rho}_{+-}`` and ``\hat{\rho}_{-+}`` between the two sectors.

In terms of Liouville spaces, the even block acts on ``\mathcal{L}(V_+ \leftarrow V_+) \oplus \mathcal{L}(V_- \leftarrow V_-)`` and the odd block on ``\mathcal{L}(V_+ \leftarrow V_-) \oplus \mathcal{L}(V_- \leftarrow V_+)``.

### Block-diagonal form of a Liouvillian

Let us add the loss of single photons (weak symmetry):

```@example symmetries
c_ops = [a, b]
L = liouvillian(H, c_ops)
bdf_L = block_diagonal_form(L)
L_even, L_odd = bdf_L.blocks
L_even
```

The sectors of the Liouvillian also tell which states of ``V`` belong together: two basis states ``|i\rangle`` and ``|i'\rangle`` belong to the same sector of ``V`` if, for every ``|j\rangle``, ``|i\rangle\langle j|`` and ``|i'\rangle\langle j|`` are in the same block of ``\mathcal{L}``, and so are ``|j\rangle\langle i|`` and ``|j\rangle\langle i'|``. Here [`block_diagonal_form`](@ref) finds ``V_+`` and ``V_-``, and describes each block of ``\mathcal{L}`` as a direct sum of Liouville spaces: `L_even` acts on `ℒ(5 ← 5) ⊕ ℒ(4 ← 4)` (``\hat{\rho}_{++}`` and ``\hat{\rho}_{--}``), and `L_odd` on `ℒ(5 ← 4) ⊕ ℒ(4 ← 5)` (``\hat{\rho}_{+-}`` and ``\hat{\rho}_{-+}``):

```@example symmetries
L_odd
```

`bdf_L.B` acts on the direct sum of all the blocks of ``\hat{\rho}``, `ℒ(5 ← 5) ⊕ ℒ(4 ← 4) ⊕ ℒ(5 ← 4) ⊕ ℒ(4 ← 5)`, since the direct sums are flattened. Which of them form a block of ``\mathcal{L}`` is given by `bdf_L.blocks`, and `bdf_L.B == direct_sum(bdf_L.blocks...)`:

```@example symmetries
bdf_L.P
```

```@example symmetries
bdf_L.B == L_even ⊕ L_odd
```

A vector of `L_even` is made of ``\hat{\rho}_{++}`` and ``\hat{\rho}_{--}``, one after the other, each vectorized as by [`mat2vec`](@ref). Therefore, [`direct_component`](@ref) and [`vec2mat`](@ref) give them back as matrices, acting on ``V_+`` (dims `([5], [5])`) and ``V_-`` (dims `([4], [4])`). For example, the steady state is the eigenvector of `L_even` with zero eigenvalue, and the trace of its block ``\hat{\rho}_{++}`` is the probability of an even number of photons:

```@example symmetries
λ, vecs = eigenstates(L_even)
v = vecs[argmin(abs.(λ))]
ρ₊₊ = vec2mat(direct_component(v, 1))
ρ₋₋ = vec2mat(direct_component(v, 2))
real(tr(ρ₊₊) / (tr(ρ₊₊) + tr(ρ₋₋)))
```

The full density matrix is obtained by bringing the vector back to the original Liouville space with `bdf_L.P'`. Since `L_even` acts on the first two summands of the space of `bdf_L.B`, `direct_sum(S, 1 => v)` puts `v` there and zeros elsewhere:

```@example symmetries
ρ = vec2mat(bdf_L.P' * direct_sum(bdf_L.B.dimensions.to, 1 => v))
ρ = ρ / tr(ρ)
ρ ≈ steadystate(H, c_ops)
```

With the loss of photons in pairs (strong symmetry), the Liouvillian has four blocks, one for each block of ``\hat{\rho}``, and there is one steady state in ``\hat{\rho}_{++}`` and one in ``\hat{\rho}_{--}``:

```@example symmetries
bdf_strong = block_diagonal_form(liouvillian(H, [a^2, b^2]))
[count(λ -> abs(λ) < 1.0e-10, eigenenergies(Lk)) for Lk in bdf_strong.blocks]
```

```@example symmetries
bdf_strong.blocks[2]
```

### Jaynes-Cummings model

With a cutoff of 3 photons, the Jaynes-Cummings Hamiltonian splits into the sectors ``N = 0, 1, 2, 3``. Here we use `destroy(2)` as ``\hat{\sigma}_-``, so that `basis(2, 0)` is the ground state ``|g\rangle``:

```@example symmetries
Nc = 3
a = destroy(Nc) ⊗ qeye(2)
σm = qeye(Nc) ⊗ destroy(2)
H_JC = a' * a + σm' * σm + 0.3 * (a' * σm + a * σm')
block_diagonal_form(H_JC).B
```

The loss of photons and the decay of the qubit both remove one excitation (weak ``U(1)`` symmetry), so they keep the difference ``d = N - N'`` between the number of excitations of the ket and of the bra of ``|\psi_N\rangle\langle\psi_{N'}|``. The Liouvillian has one block for each ``d = -3, \ldots, 3``, and the block ``d`` is made of the blocks ``\hat{\rho}_{N, N'}`` with ``N - N' = d``:

```@example symmetries
bdf_JC = block_diagonal_form(liouvillian(H_JC, [a, σm]))
bdf_JC.block_sizes
```

```@example symmetries
bdf_JC.blocks[1] # d = 0: the populations of each sector and the coherences inside it
```

```@example symmetries
bdf_JC.blocks[2] # d = 1
```

## Summary

| System | Blocks of the Hamiltonian | Blocks of the Liouvillian |
|:-------|:--------------------------|:--------------------------|
| Two modes, weak ``\mathbb{Z}_2`` (``\hat{a}``, ``\hat{b}``) | ``V_+``, ``V_-`` | ``\mathcal{L}(V_+ \leftarrow V_+) \oplus \mathcal{L}(V_- \leftarrow V_-)``, ``\mathcal{L}(V_+ \leftarrow V_-) \oplus \mathcal{L}(V_- \leftarrow V_+)`` |
| Two modes, strong ``\mathbb{Z}_2`` (``\hat{a}^2``, ``\hat{b}^2``) | ``V_+``, ``V_-`` | ``\mathcal{L}(V_+ \leftarrow V_+)``, ``\mathcal{L}(V_- \leftarrow V_+)``, ``\mathcal{L}(V_+ \leftarrow V_-)``, ``\mathcal{L}(V_- \leftarrow V_-)`` |
| Jaynes-Cummings, weak ``U(1)`` (``\hat{a}``, ``\hat{\sigma}_-``) | ``V_0``, ``V_1``, … | ``\bigoplus_{N - N' = d} \mathcal{L}(V_N \leftarrow V_{N'})`` for each ``d`` |

!!! note
    [`block_diagonal_form`](@ref) only finds the blocks of a matrix in the basis it is written in. It finds the symmetries that are diagonal in that basis (such as the parity or the number of excitations in the Fock basis), but not, e.g., the exchange of two identical modes, whose sectors are the symmetric and antisymmetric combinations of Fock states.
