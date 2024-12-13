```@meta
CurrentModule = QuantumToolbox
```

# [API](@id doc-API)

**Table of contents**

[[toc]] <!-- the level setting is in ".vitepress/config.mts" -->

## [Quantum object (Qobj) and type](@id doc-API:Quantum-object-and-type)

```@docs
AbstractQuantumObject
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
QuantumObjectEvolution
Base.size
Base.eltype
Base.length
SciMLOperators.cache_operator
```

## [Qobj boolean functions](@id doc-API:Qobj-boolean-functions)

```@docs
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
SciMLOperators.iscached
SciMLOperators.isconstant
```

## [Qobj arithmetic and attributes](@id doc-API:Qobj-arithmetic-and-attributes)

```@docs
Base.zero
Base.one
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
QobjEvo
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
TimeEvolutionProblem
TimeEvolutionSol
TimeEvolutionMCSol
TimeEvolutionSSESol
sesolveProblem
mesolveProblem
mcsolveProblem
mcsolveEnsembleProblem
ssesolveProblem
ssesolveEnsembleProblem
sesolve
mesolve
mcsolve
ssesolve
dfd_mesolve
liouvillian
liouvillian_generalized
```

### [Steady State Solvers](@id doc-API:Steady-State-Solvers)

```@docs
steadystate
steadystate_floquet
SteadyStateDirectSolver
SteadyStateEigenSolver
SteadyStateLinearSolver
SteadyStateODESolver
```

### [Dynamical Shifted Fock method](@id doc-API:Dynamical-Shifted-Fock-method)

```@docs
dsf_mesolve
dsf_mcsolve
```

### [Low-rank time evolution](@id doc-API:Low-rank-time-evolution)

```@docs
TimeEvolutionLRSol
lr_mesolveProblem
lr_mesolve
```

## [Correlations and Spectrum](@id doc-API:Correlations-and-Spectrum)

```@docs
correlation_3op_2t
correlation_3op_1t
correlation_2op_2t
correlation_2op_1t
spectrum_correlation_fft
spectrum
ExponentialSeries
PseudoInverse
```

## [Metrics](@id doc-API:Metrics)

```@docs
entropy_vn
entanglement
tracedist
fidelity
```

## [Spin Lattice](@id doc-API:Spin-Lattice)

```@docs
Lattice
MultiSiteOperator
DissipativeIsing
```

## [Miscellaneous](@id doc-API:Miscellaneous)

```@docs
wigner
negativity
```

## [Linear Maps](@id doc-API:Linear-Maps)

```@docs
AbstractLinearMap
```

## [Utility functions](@id doc-API:Utility-functions)

```@docs
QuantumToolbox.versioninfo
QuantumToolbox.about
gaussian
n_thermal
PhysicalConstants
convert_unit
row_major_reshape
meshgrid
```

## [Visualization](@id doc-API:Visualization)

```@docs
plot_wigner
```
