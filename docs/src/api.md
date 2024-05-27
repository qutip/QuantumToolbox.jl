```@meta
CurrentModule = QuantumToolbox
```

# [API](@id doc-API)

## Contents

```@contents
Pages = ["api.md"]
```

## [Quantum object functions](@id doc-API: Quantum object functions)

```@docs
BraQuantumObject
Bra
KetQuantumObject
Ket
OperatorQuantumObject
Operator
OperatorBraQuantumObject
OperatorBra
OperatorKetQuantumObject
OperatorKet
SuperOperatorQuantumObject
SuperOperator
QuantumObject
Qobj
ket2dm
isbra
isket
isoper
isoperbra
isoperket
issuper
size
eltype
length
sqrtm
LinearAlgebra.tr
LinearAlgebra.svdvals
LinearAlgebra.norm
LinearAlgebra.kron
tensor
⊗
```

## [General functions](@id doc-API: General functions)

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
fidelity
wigner
get_coherence
n_th
tracedist
get_data
mat2vec
vec2mat
```

## [Quantum states, operators and super-operators](@id doc-API: Quantum states, operators and super-operators)

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
qeye
fock
basis
coherent
rand_dm
projection
sinm
cosm
```

## [Time evolution](@id doc-API: Time evolution)

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

## [Correlations and Spectrum](@id doc-API: Correlations and Spectrum)

```@docs
correlation_3op_2t
correlation_2op_2t
correlation_2op_1t
spectrum
```

## [Eigenvalues and eigenvectors](@id doc-API: Eigenvalues and eigenvectors)

```@docs
EigsolveResult
eigenenergies
eigenstates
LinearAlgebra.eigen
LinearAlgebra.eigvals
eigsolve
eigsolve_al
```

## [Low Rank internal APIs](@id doc-API: Low Rank internal APIs)

```@docs
_calculate_expectation!
_adjM_condition_variational
_adjM_affect!
_adjM_condition_ratio
_pinv!
dBdz!
```

## [Miscellaneous](@id doc-API: Miscellaneous)

```@docs
QuantumToolbox.versioninfo
QuantumToolbox.about
```
