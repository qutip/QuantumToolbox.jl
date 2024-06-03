```@meta
CurrentModule = QuantumToolbox
```

# [API](@id doc-API)

## Contents

```@contents
Pages = ["api.md"]
```

## [Quantum object (Qobj) and type](@id doc-API:Quantum-object-and-type)

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
size
eltype
length
isbra
isket
isoper
isoperbra
isoperket
issuper
LinearAlgebra.ishermitian
isherm
LinearAlgebra.issymmetric
LinearAlgebra.isposdef
```

## [Qobj arithmetic and attributes](@id doc-API:Qobj-arithmetic-and-attributes)

```@docs
Base.conj
LinearAlgebra.transpose
trans
LinearAlgebra.adjoint
dag
dagger
matrix_element
LinearAlgebra.sqrt
sqrtm
sinm
cosm
LinearAlgebra.tr
LinearAlgebra.svdvals
LinearAlgebra.norm
LinearAlgebra.normalize
LinearAlgebra.normalize!
unit
LinearAlgebra.inv
ptrace
tidyup
tidyup!
get_data
get_coherence
partial_transpose
```

## [Qobj eigenvalues and eigenvectors](@id doc-API:Qobj-eigenvalues-and-eigenvectors)

```@docs
EigsolveResult
eigenenergies
eigenstates
LinearAlgebra.eigen
LinearAlgebra.eigvals
eigsolve
eigsolve_al
```

## [Qobj manipulation](@id doc-API:Qobj-manipulation)

```@docs
ket2dm
expect
LinearAlgebra.kron
tensor
⊗
sparse_to_dense
dense_to_sparse
vec2mat
mat2vec
```

## [Generate states and operators](@id doc-API:Generate-states-and-operators)

```@docs
fock
basis
coherent
rand_dm
sigmap
sigmam
sigmax
sigmay
sigmaz
destroy
create
eye
qeye
projection
spre
spost
sprepost
lindblad_dissipator
```

## [Time evolution](@id doc-API:Time-evolution)

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

## [Correlations and Spectrum](@id doc-API:Correlations-and-Spectrum)

```@docs
correlation_3op_2t
correlation_2op_2t
correlation_2op_1t
spectrum
```

## [Metrics](@id doc-API:Metrics)

```@docs
entropy_vn
entanglement
tracedist
fidelity
```

## [Miscellaneous](@id doc-API:Miscellaneous)

```@docs
wigner
negativity
```

## [Utility functions](@id doc-API:Utility-functions)

```@docs
QuantumToolbox.versioninfo
QuantumToolbox.about
gaussian
n_th
row_major_reshape
meshgrid
_calculate_expectation!
_adjM_condition_variational
_adjM_affect!
_adjM_condition_ratio
_pinv!
dBdz!
```
