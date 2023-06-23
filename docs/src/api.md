```@meta
CurrentModule = QuPhys
```

# [API](@id API)

## [Quantum object functions](@id API: Quantum object functions)

```@docs
BraQuantumObject
KetQuantumObject
OperatorQuantumObject
SuperOperatorQuantumObject
QuantumObject
ket2dm
isbra
isket
isoper
issuper
size
length
LinearAlgebra.tr
LinearAlgebra.norm
LinearAlgebra.kron
LinearAlgebra.eigen
LinearAlgebra.eigvals
```

## [General functions](@id API: General functions)

```@docs
row_major_reshape
meshgrid
sparse_to_dense
dense_to_sparse
gaussian
ptrace
entropy_vn
entanglement
expect
wigner
get_coherence
n_th
get_data
mat2vec
vec2mat
```

## [Quantum states, operators and super-operators](@id API: Quantum states, operators and super-operators)

```@docs
spre
spost
sprepost
lindblad_dissipator
destroy
create
sigmap
sigmam
sigmax
sigmay
sigmaz
eye
fock
basis
coherent
rand_dm
projection
sinm
cosm
```

## [Time evolution](@id API: Time evolution)

```@docs
sesolveProblem
mesolveProblem
mcsolveProblem
mcsolveEnsembleProblem
sesolve
mesolve
mcsolve
dfd_mesolve
dsf_mesolve
dsf_mcsolve
liouvillian
liouvillian_generalized
steadystate_floquet
```

## [Correlations and Spectrum](@id API: Correlations and Spectrum)
```@docs
correlation_3op_2t
correlation_2op_2t
correlation_2op_1t
spectrum
```