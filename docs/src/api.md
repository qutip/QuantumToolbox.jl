```@meta
CurrentModule = QuantumToolbox
```

# [API](@id API)

## [Quantum object functions](@id API: Quantum object functions)

```@docs
BraQuantumObject
KetQuantumObject
OperatorQuantumObject
SuperOperatorQuantumObject
QuantumObject
Qobj
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
tensor
⊗
```

## [General functions](@id API: General functions)

```@docs
row_major_reshape
meshgrid
sparse_to_dense
dense_to_sparse
tidyup
tidyup!
gaussian
ptrace
partial_transpose
negativity
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
lr_mesolveProblem
mcsolveProblem
mcsolveEnsembleProblem
sesolve
mesolve
lr_mesolve
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

## [Eigenvalues and eigenvectors](@id API: Eigenvalues and eigenvectors)
```@docs
LinearAlgebra.eigen
LinearAlgebra.eigvals
eigsolve
eigsolve_al
```

## [Low Rank internal APIs](@id API: Low Rank internal APIs)
```@docs
_calculate_expectation!
_adjM_condition_variational
_adjM_affect!
_adjM_condition_ratio
_pinv!
dBdz!
```
