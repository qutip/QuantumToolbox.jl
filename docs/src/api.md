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
OperatorSum
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
LinearAlgebra.issymmetric
LinearAlgebra.isposdef
isunitary
```

## [Qobj arithmetic and attributes](@id doc-API:Qobj-arithmetic-and-attributes)

```@docs
Base.conj
LinearAlgebra.transpose
LinearAlgebra.adjoint
LinearAlgebra.dot
LinearAlgebra.sqrt
LinearAlgebra.log
LinearAlgebra.exp
LinearAlgebra.sin
LinearAlgebra.cos
LinearAlgebra.tr
LinearAlgebra.svdvals
LinearAlgebra.norm
LinearAlgebra.normalize
LinearAlgebra.normalize!
LinearAlgebra.inv
LinearAlgebra.diag
proj
ptrace
purity
permute
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
variance
LinearAlgebra.kron
sparse_to_dense
dense_to_sparse
vec2mat
mat2vec
```

## [Generate states and operators](@id doc-API:Generate-states-and-operators)

```@docs
zero_ket
fock
basis
coherent
rand_ket
fock_dm
coherent_dm
thermal_dm
maximally_mixed_dm
rand_dm
spin_state
spin_coherent
bell_state
singlet_state
triplet_states
w_state
ghz_state
rand_unitary
sigmap
sigmam
sigmax
sigmay
sigmaz
jmat
spin_Jx
spin_Jy
spin_Jz
spin_Jm
spin_Jp
spin_J_set
destroy
create
displace
squeeze
num
QuantumToolbox.position
QuantumToolbox.momentum
phase
fdestroy
fcreate
tunneling
qft
eye
projection
commutator
spre
spost
sprepost
lindblad_dissipator
```

## [Synonyms of functions for Qobj](@id doc-API:Synonyms-of-functions-for-Qobj)
```@docs
Qobj
shape
isherm
trans
dag
matrix_element
unit
sqrtm
logm
expm
sinm
cosm
tensor
⊗
qeye
```

## [Time evolution](@id doc-API:Time-evolution)

```@docs
TimeEvolutionSol
TimeEvolutionMCSol
sesolveProblem
mesolveProblem
lr_mesolveProblem
mcsolveProblem
ssesolveProblem
mcsolveEnsembleProblem
ssesolveEnsembleProblem
sesolve
mesolve
lr_mesolve
mcsolve
ssesolve
dfd_mesolve
dsf_mesolve
dsf_mcsolve
liouvillian
liouvillian_generalized
steadystate
steadystate_floquet
SteadyStateODESolver
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
