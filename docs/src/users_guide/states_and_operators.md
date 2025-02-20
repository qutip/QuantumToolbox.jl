# [Manipulating States and Operators](@id doc:Manipulating-States-and-Operators)

## Introduction
In the previous guide section [Basic Operations on Quantum Objects](@ref doc:Qobj), we saw how to create states and operators, using the functions built into `QuantumToolbox`. In this portion of the guide, we will look at performing basic operations with states and operators. For more detailed demonstrations on how to use and manipulate these objects, see the examples given in the tutorial section.

```@setup states_and_operators
using QuantumToolbox
```

## [State Vectors (kets or bras)](@id doc:State-vectors)
Here we begin by creating a Fock [`basis`](@ref) (or [`fock`](@ref)) vacuum state vector ``|0\rangle`` with in a Hilbert space with `5` number states, from `0` to `4`:

```@example states_and_operators
vac = basis(5, 0)
```

and then create a lowering operator ``\hat{a}`` corresponding to `5` number states using the [`destroy`](@ref) function:

```@example states_and_operators
a = destroy(5)
```

Now lets apply the lowering operator ``\hat{a}`` to our vacuum state `vac`:

```@example states_and_operators
a * vac
```

We see that, as expected, the vacuum is transformed to the zero vector. A more interesting example comes from using the `adjoint` of the lowering operator ``\hat{a}``, the raising operator ``\hat{a}^\dagger``:

```@example states_and_operators
a' * vac
```

The raising operator has in indeed raised the state `vac` from the vacuum to the ``|1\rangle`` state. Instead of using the `adjoint` method to raise the state, we could have also used the built-in [`create`](@ref) function to make a raising operator:

```@example states_and_operators
ad = create(5)
ad * vac
```

which does the same thing. We can raise the vacuum state more than once by successively apply the raising operator:

```@example states_and_operators
ad * ad * vac
```

or just taking the square of the raising operator ``\left(\hat{a}^\dagger\right)^2``:

```@example states_and_operators
ad^2 * vac
```

Applying the raising operator twice gives the expected ``\sqrt{n+1}`` dependence. We can use the product of ``\hat{a}^\dagger \hat{a}`` to also apply the number operator to the state vector `vac`:

```@example states_and_operators
ad * a * vac
```

or on the ``|1\rangle`` state:

```@example states_and_operators
ad * a * (ad * vac)
```

or on the ``|2\rangle`` state:

```@example states_and_operators
ad * a * (ad^2 * vac)
```

Notice how in this last example, application of the number operator does not give the expected value ``n=2``, but rather ``2\sqrt{2}``. This is because this last state is not normalized to unity as ``\hat{a}^\dagger|n\rangle=\sqrt{n+1}|n+1\rangle``. Therefore, we should [`normalize`](@ref) (or use [`unit`](@ref)) our vector first:

```@example states_and_operators
ad * a * normalize(ad^2 * vac)
```

Since we are giving a demonstration of using states and operators, we have done a lot more work than we should have. For example, we do not need to operate on the vacuum state to generate a higher number Fock state. Instead we can use the [`basis`](@ref) (or [`fock`](@ref)) function to directly obtain the required state:

```@example states_and_operators
ket = basis(5, 2)
```

Notice how it is automatically normalized. We can also use the built in number operator [`num`](@ref):

```@example states_and_operators
n = num(5)
```

Therefore, instead of `ad * a * normalize(ad^2 * vac)`, we have:

```@example states_and_operators
n * ket
```

We can also create superpositions of states:

```@example states_and_operators
ket = normalize(basis(5, 0) + basis(5, 1))
```

where we have used the `normalize` function again to normalize the state. Apply the number operator again:

```@example states_and_operators
n * ket
```

We can also create coherent states and squeezed states by applying the [`displace`](@ref) and [`squeeze`](@ref) functions to the vacuum state:

```@example states_and_operators
vac = basis(5, 0)

d = displace(5, 1im)

s = squeeze(5, 0.25 + 0.25im)

d * vac
```

```@example states_and_operators
d * s * vac
```

Of course, displacing the vacuum gives a coherent state, which can also be generated using the built in [`coherent`](@ref) function.

## [Density matrices](@id doc:Density-matrices)

One of the main purpose of `QuantumToolbox` is to explore the dynamics of open quantum systems, where the most general state of a system is no longer a state vector, but rather a density matrix. Since operations on density matrices operate identically to those of vectors, we will just briefly highlight creating and using these structures.

The simplest density matrix is created by forming the outer-product ``|\psi\rangle\langle\psi|`` of a ket vector:

```@example states_and_operators
ket = basis(5, 2)
ket * ket'
```

A similar task can also be accomplished via the [`fock_dm`](@ref) or [`ket2dm`](@ref) functions:

```@example states_and_operators
fock_dm(5, 2)
```

```@example states_and_operators
ket2dm(ket)
```

If we want to create a density matrix with equal classical probability of being found in the ``|2\rangle`` or ``|4\rangle`` number states, we can do the following:

```@example states_and_operators
0.5 * fock_dm(5, 2) + 0.5 * fock_dm(5, 4) # with fock_dm
0.5 * ket2dm(basis(5, 2)) + 0.5 * ket2dm(basis(5, 4)) # with ket2dm
```

There are also several other built-in functions for creating predefined density matrices, for example [`coherent_dm`](@ref) and [`thermal_dm`](@ref) which create coherent state and thermal state density matrices, respectively.

```@example states_and_operators
coherent_dm(5, 1.25)
```

```@example states_and_operators
thermal_dm(5, 1.25)
```

`QuantumToolbox` also provides a set of distance metrics (see section [Entropy and Metrics](@ref doc-API:Entropy-and-Metrics) in API page) for determining how close two density matrix distributions are to each other. Included are the [`fidelity`](@ref), and trace distance ([`tracedist`](@ref)).

```@example states_and_operators
x = coherent_dm(5, 1.25)

y = coherent_dm(5, 1.25im)

z = thermal_dm(5, 0.125)

fidelity(x, y)
```
Note that the definition of [`fidelity`](@ref) here is from [Nielsen-Chuang2011](@citet). It is the square root of the fidelity defined in [Jozsa1994](@citet). We also know that for two pure states, the trace distance (``T``) and the fidelity (``F``) are related by ``T = \sqrt{1-F^2}``:

```@example states_and_operators
tracedist(x, y) ≈ sqrt(1 - (fidelity(x, y))^2)
```

For a pure state and a mixed state, ``1 - F \leq T`` which can also be verified:

```@example states_and_operators
1 - fidelity(x, z) < tracedist(x, z)
```

## [Two-level systems (Qubits)](@id doc:Two-level-systems)

Having spent a fair amount of time on basis states that represent harmonic oscillator states, we now move on to qubit, or two-level quantum systems (for example a spin-``1/2``). To create a state vector corresponding to a qubit system, we use the same basis, or fock, function with only two levels:

```@example states_and_operators
spin = basis(2, 0)
```

Now at this point one may ask how this state is different than that of a harmonic oscillator in the vacuum state truncated to two energy levels?

```@example states_and_operators
vac = basis(2, 0)
```

At this stage, there is no difference. This should not be surprising as we called the exact same function twice. The difference between the two comes from the action of the spin operators [`sigmax`](@ref), [`sigmay`](@ref), [`sigmaz`](@ref), [`sigmap`](@ref), and [`sigmam`](@ref) on these two-level states. For example, if `vac` corresponds to the vacuum state of a harmonic oscillator, then, as we have already seen, we can use the raising operator ([`create`](@ref)) to get the ``|1\rangle`` state:

```@example states_and_operators
create(2) * vac
```

For a spin system, the operator analogous to the raising operator is the ``\hat{\sigma}_+`` operator [`sigmap`](@ref). Applying on the spin state gives:

```@example states_and_operators
sigmap() * spin
```

Now we see the difference! The [`sigmap`](@ref) operator acting on the spin state returns the zero vector. Why is this? To see what happened, let us use the ``\hat{\sigma}_z`` ([`sigmaz`](@ref)) operator:

```@example states_and_operators
sigmaz()
```

```@example states_and_operators
sigmaz() * spin
```

```@example states_and_operators
spin2 = basis(2, 1)
```

```@example states_and_operators
sigmaz() * spin2
```

The answer is now apparent. Since the `QuantumToolbox` [`sigmaz`](@ref) function uses the standard ``Z``-basis representation of the ``\hat{\sigma}_z`` spin operator, the `spin` state corresponds to the ``|\uparrow\rangle`` state of a two-level spin system while `spin2` gives the ``|\downarrow\rangle`` state. Therefore, in our previous example `sigmap() * spin`, we raised the qubit state out of the truncated two-level Hilbert space resulting in the zero state.

While at first glance this convention might seem somewhat odd, it is in fact quite handy. For one, the spin operators remain in the conventional form. Second, this corresponds nicely with the quantum information definitions of qubit states, where the excited ``|\uparrow\rangle`` state is label as ``|0\rangle``, and the ``|\downarrow\rangle`` state by ``|1\rangle``.

If one wants to create spin operators for higher spin systems, then the [`jmat`](@ref) function comes in handy.

## [Expectation values](@id doc:Expectation-values)

Some of the most important information about quantum systems comes from calculating the expectation value of operators, both Hermitian and non-Hermitian, as the state or density matrix of the system varies in time. Therefore, in this section we demonstrate the use of the [`expect`](@ref) function. To begin:

```@example states_and_operators
vac = basis(5, 0)

one = basis(5, 1)

c = create(5)

N = num(5)

coh = coherent_dm(5, 1.0im)

cat = normalize(basis(5, 4) + 1.0im * basis(5, 3))

println(expect(N, vac) ≈ 0)
println(expect(N, one) ≈ 1)
println(expect(N, coh) ≈ 0.9970555745806597)
println(expect(c, cat) ≈ 1im)
```

The [`expect`](@ref) function also accepts lists or arrays of state vectors or density matrices for the second input:

```@example states_and_operators
states = [normalize(c^k * vac) for k in 0:4]

expect(N, states)
```

```@example states_and_operators
cat_list = [normalize(basis(5, 4) + x * basis(5, 3)) for x in [0, 1.0im, -1.0, -1.0im]]

expect(c, cat_list)
```

Notice how in this last example, all of the return values are complex numbers. This is because the expect function looks to see whether the operator is Hermitian or not. If the operator is Hermitian, then the output will always be real. In the case of non-Hermitian operators, the return values may be complex. Therefore, the expect function will return an array of complex values for non-Hermitian operators when the input is a list/array of states or density matrices.

Of course, the expect function works for spin states and operators:

```@example states_and_operators
up = basis(2, 0)

dn = basis(2, 1)

println(expect(sigmaz(), up) ≈ 1)
println(expect(sigmaz(), dn) ≈ -1)
```

as well as the composite objects discussed in the next section [Tensor Products and Partial Traces](@ref doc:Tensor-products-and-Partial-Traces):

```@example states_and_operators
spin1 = basis(2, 0)

spin2 = basis(2, 1)

two_spins = tensor(spin1, spin2)

sz1 = tensor(sigmaz(), qeye(2))

sz2 = tensor(qeye(2), sigmaz())

println(expect(sz1, two_spins) ≈ 1)
println(expect(sz2, two_spins) ≈ -1)
```

## [Superoperators and Vectorized Operators](@id doc:Superoperators-and-Vectorized-Operators)

In addition to state vectors and density operators, `QuantumToolbox` allows for representing maps that act linearly on density operators using the Liouville supermatrix formalisms.

This support is based on the correspondence between linear operators acting on a Hilbert space, and vectors in two copies of that Hilbert space (which is also called the Fock-Liouville space), 
```math
\textrm{vec} : \mathcal{L}(\mathcal{H}) \rightarrow \mathcal{H}\otimes\mathcal{H}.
```
Therefore, a given density matrix ``\hat{\rho}`` can then be vectorized, denoted as 
```math
|\hat{\rho}\rangle\rangle = \textrm{vec}(\hat{\rho}).
```

`QuantumToolbox` uses the column-stacking convention for the isomorphism between ``\mathcal{L}(\mathcal{H})`` and ``\mathcal{H}\otimes\mathcal{H}``. This isomorphism is implemented by the functions [`mat2vec`](@ref) (or [`operator_to_vector`](@ref)) and [`vec2mat`](@ref) (or [`vector_to_operator`](@ref)):

```@example states_and_operators
rho = Qobj([1 2; 3 4])
```

```@example states_and_operators
vec_rho = mat2vec(rho)
```

```@example states_and_operators
rho2 = vec2mat(vec_rho)
```

The `QuantumObject.type` attribute indicates whether a quantum object is a vector corresponding to an [`OperatorKet`](@ref), or its Hermitian conjugate [`OperatorBra`](@ref). One can also use [`isoper`](@ref), [`isoperket`](@ref), and [`isoperbra`](@ref) to check the type:

```@example states_and_operators
println(isoper(vec_rho))
println(isoperket(vec_rho))
println(isoperbra(vec_rho))
println(isoper(vec_rho'))
println(isoperket(vec_rho'))
println(isoperbra(vec_rho'))
```

Because `Julia` is a column-oriented languages (like `Fortran` and `MATLAB`), in `QuantumToolbox`, we define the [`spre`](@ref) (left), [`spost`](@ref) (right), and [`sprepost`](@ref) (left-and-right) multiplication superoperators as follows:

```math
\begin{align}
\hat{A}\hat{\rho}~~~ &\rightarrow \textrm{spre}(\hat{A}) * \textrm{vec}(\hat{\rho}) = \hat{\mathbb{1}}\otimes \hat{A} ~ |\hat{\rho}\rangle\rangle,\notag\\
\hat{\rho} \hat{B} &\rightarrow \textrm{spost}(\hat{B}) * \textrm{vec}(\hat{\rho}) = \hat{B}^T\otimes \hat{\mathbb{1}} ~ |\hat{\rho}\rangle\rangle,\notag\\
\hat{A} \hat{\rho} \hat{B} &\rightarrow \textrm{sprepost}(\hat{A},\hat{B}) * \textrm{vec}(\hat{\rho}) = \hat{B}^T\otimes \hat{A} ~ |\hat{\rho}\rangle\rangle,\notag
\end{align}
```
where ``\hat{\mathbb{1}}`` represents the identity operator with Hilbert space dimension equal to ``\hat{\rho}``.

```@example states_and_operators
A = Qobj([1 2; 3 4])
S_A = spre(A)
```

```@example states_and_operators
B = Qobj([5 6; 7 8])
S_B = spost(B)
```

```@example states_and_operators
S_AB = sprepost(A, B)
```

```@example states_and_operators
S_AB ≈ S_A * S_B ≈ S_B * S_A
```

One can also use [`issuper`](@ref) to check the type:

```@example states_and_operators
println(isoper(S_AB))
println(issuper(S_AB))
```

With the above definitions, the following equalities hold in `Julia`:

```math
\textrm{vec}(\hat{A} \hat{\rho} \hat{B}) = \textrm{spre}(\hat{A}) * \textrm{spre}(\hat{B}) * \textrm{vec}(\hat{\rho}) = \textrm{sprepost}(\hat{A},\hat{B}) * \textrm{vec}(\hat{\rho}) ~~\forall~~\hat{A}, \hat{B}, \hat{\rho}
```

```@example states_and_operators
N  = 10
A = Qobj(rand(ComplexF64, N, N))
B = Qobj(rand(ComplexF64, N, N))
ρ = rand_dm(N) # random density matrix
mat2vec(A * ρ * B) ≈ spre(A) * spost(B) * mat2vec(ρ) ≈ sprepost(A, B) * mat2vec(ρ)
```

In addition, dynamical generators on this extended space, often called Liouvillian superoperators, can be created using the [`liouvillian`](@ref) function. Each of these takes a Hamiltonian along with a list of collapse operators, and returns a [`type=SuperOperator`](@ref SuperOperator) object that can be exponentiated to find the superoperator for that evolution.

```@example states_and_operators
H = 10 * sigmaz()

c = destroy(2)

L = liouvillian(H, [c])
```

```@example states_and_operators
t = 0.8
exp(L * t)
```

See the section [Lindblad Master Equation Solver](@ref doc-TE:Lindblad-Master-Equation-Solver) for more details.
