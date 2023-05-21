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
gaussian
ptrace
entropy_vn
entanglement
expect
wigner
get_coherence
n_th
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
liouvillian
liouvillian_generalized
sesolve
mesolve
mcsolve
dfd_mesolve
dsf_mesolve
dsf_mcsolve
```

## [Correlations and Spectrum](@id API: Correlations and Spectrum)
```@docs
correlation_3op_2t
correlation_2op_2t
correlation_2op_1t
spectrum
```