```@meta
CurrentModule = QuantumToolbox

DocTestSetup = quote
    using LinearAlgebra
    using SparseArrays
    using QuantumToolbox
end
```

# [API](@id doc-API)

**Table of contents**

[[toc]] <!-- the level setting is in ".vitepress/config.mts" -->

## [Quantum object (Qobj) and type](@id doc-API:Quantum-object-and-type)

```@docs
AbstractSpace
HilbertSpace
EnrSpace
ProductDimensions
AbstractQuantumObject
Bra
Ket
Operator
OperatorBra
OperatorKet
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
SparseArrays.permute
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
to_dense
to_sparse
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
enr_fock
enr_thermal_dm
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
enr_destroy
enr_identity
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
tensor
⊗
qeye
vector_to_operator
operator_to_vector
sqrtm
logm
expm
sinm
cosm
qeye_like
qzero_like
```

## [Time evolution](@id doc-API:Time-evolution)

```@docs
TimeEvolutionProblem
TimeEvolutionSol
TimeEvolutionMCSol
TimeEvolutionStochasticSol
average_states
average_expect
std_expect
sesolveProblem
mesolveProblem
mcsolveProblem
mcsolveEnsembleProblem
ssesolveProblem
ssesolveEnsembleProblem
smesolveProblem
smesolveEnsembleProblem
sesolve
mesolve
mcsolve
ssesolve
smesolve
sesolve_map
mesolve_map
dfd_mesolve
liouvillian
liouvillian_dressed_nonsecular
bloch_redfield_tensor
brterm
brmesolve
```

### [Steady State Solvers](@id doc-API:Steady-State-Solvers)

```@docs
steadystate
steadystate_fourier
SteadyStateDirectSolver
SteadyStateEigenSolver
SteadyStateLinearSolver
SteadyStateODESolver
SSFloquetEffectiveLiouvillian
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
Lanczos
```

## [Entropy and Metrics](@id doc-API:Entropy-and-Metrics)

```@docs
entropy_vn
entropy_relative
entropy_linear
entropy_mutual
entropy_conditional
entanglement
concurrence
negativity
fidelity
tracedist
hilbert_dist
hellinger_dist
bures_dist
bures_angle
```

## [Spin Lattice](@id doc-API:Spin-Lattice)

```@docs
Lattice
multisite_operator
DissipativeIsing
```

## [Symmetries and Block Diagonalization](@id doc-API:Symmetries-and-Block-Diagonalization)

```@docs
block_diagonal_form
BlockDiagonalForm
```

## [Miscellaneous](@id doc-API:Miscellaneous)

```@docs
wigner
```

## [Linear Maps](@id doc-API:Linear-Maps)

```@docs
AbstractLinearMap
```

## [Utility functions](@id doc-API:Utility-functions)

```@docs
QuantumToolbox.settings
QuantumToolbox.versioninfo
QuantumToolbox.about
QuantumToolbox.cite
gaussian
n_thermal
PhysicalConstants
convert_unit
row_major_reshape
meshgrid
enr_state_dictionaries
get_hilbert_size
get_liouville_size
```

## [Visualization](@id doc-API:Visualization)

```@docs
plot_wigner
plot_fock_distribution
matrix_heatmap
matrix_histogram
```

### [Bloch Sphere](@id doc-API:Bloch-Sphere)

```@docs
Bloch
plot_bloch
render
add_points!
add_vectors!
add_line!
add_arc!
add_states!
clear!
```
