# ChangeLog

All notable changes to [`QuantumToolbox.jl`](https://github.com/qutip/QuantumToolbox.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/qutip/QuantumToolbox.jl/tree/main)

## [v0.35.0]
Release date: 2025-09-03

- Add support of `QobjEvo` for `steadystate` (ODE solver only). ([#536])
- Changes to `SteadyStateODESolver`. ([#537])
  - Introduce the tolerances for `steadystate` terminate condition (two new fields: `terminate_reltol = 1e-5` and `terminate_abstol = 1e-7`)
  - Fix keyword argument handling for `SteadyStateODESolver` before passing to `mesolve`.
- Fix incorrect `negativity` and `partial_transpose` for arbitrary subsystem dimension. ([#539])

## [v0.34.1]
Release date: 2025-08-23

- Improve Bloch sphere rendering for animation. ([#520])
- Add support to `Enzyme.jl` for `sesolve` and `mesolve`. ([#531])

## [v0.34.0]
Release date: 2025-07-29

- Improve efficiency of `bloch_redfield_tensor` by avoiding unnecessary conversions. ([#509])
- Support `SciMLOperators v1.4+`. ([#470])
- Fix compatibility with `Makie v0.24+`. ([#513])
- Add `keep_runs_results` option for multi-trajectory solvers to align with `QuTiP`. ([#512])
  - Breaking changes for multi-trajectory solutions:
    - the original fields `expect` and `states` now store the results depend on keyword argument `keep_runs_results` (decide whether to store all trajectories results or not).
    - remove field `average_expect`
    - remove field `runs_expect`
  - New statistical analysis functions:
    - `average_states`
    - `average_expect`
    - `std_expect`
- Add support to ForwardDiff.jl for `sesolve` and `mesolve`. ([#515])
- Add documentation about automatic differentiation. ([#517])

## [v0.33.0]
Release date: 2025-07-22

- Implement `EnrSpace` and corresponding functionality. ([#500])
- Check for orthogonality breakdown in `Lanczos` solver for `spectrum`. ([#501])
- Store both `times` and `times_states` in time evolution solutions. ([#506], [#504])
- Fix errors in `Julia v1.12`. ([#507])

## [v0.32.1]
Release date: 2025-06-24

This is a release just for updating documentation.

## [v0.32.0]
Release date: 2025-06-23

- Introduce `Lanczos` solver for `spectrum`. ([#476])
- Add Bloch-Redfield master equation solver. ([#473])
- Implement Bloch Sphere rendering and align style with qutip. ([#472], [#480], [#485], [#487], [#489])
- Add `Base.copy` method for `AbstractQuantumObject`. ([#486])
- Add documentation for Bloch-Redfield master equation. ([#494])

## [v0.31.1]
Release date: 2025-05-16

- Introduce `QuantumToolbox.settings` and `auto_tidyup`. ([#460])

## [v0.31.0]
Release date: 2025-05-03

- Return `sesolve` when `mesolve` allows it. ([#455])
- Simplify structure of `QuantumObjectType`s. ([#456])

## [v0.30.1]
Release date: 2025-04-24

- Support different length for `to` and `from` on GeneralDimensions. ([#448])
- Extend the `Makie.jl` extension to all the other available backends. ([#450])
- Fix definition of noise derivative in stochastic solvers. ([#453])

## [v0.30.0]
Release date: 2025-04-12

- Make CUDA conversion more general using Adapt.jl. ([#436], [#437])
- Make the generation of `fock` states non-mutating to support Zygote.jl. ([#438])
- Remove Reexport.jl from the dependencies. ([#443])
- Add support for automatic differentiation for `sesolve` and `mesolve`. ([#440])

## [v0.29.1]
Release date: 2025-03-07

- Minor changes for GPU matrices element type and word size handling. ([#430])

## [v0.29.0]
Release date: 2025-03-07

- Add support for `OperatorKet` state input for `mesolve` and `smesolve`. ([#423])
- Introduce `plot_fock_distribution` to plot the population of a state (ket, bra, or density matrix) in its basis (assumed to be Fock basis). ([#428])

## [v0.28.0]
Release date: 2025-02-22

- Support for single `AbstractQuantumObject` in `sc_ops` for faster specific method in `ssesolve` and `smesolve`. ([#408])
- Change save callbacks from `PresetTimeCallback` to `FunctionCallingCallback`. ([#410])
- Align `eigenstates` and `eigenenergies` to QuTiP. ([#411])
- Introduce `vector_to_operator` and `operator_to_vector`. ([#413])
- Introduce some entropy related functions. ([#414], [#416])
  - `entropy_linear`
  - `entropy_mutual`
  - `entropy_conditional`
  - `entropy_relative`
- Fix `entanglement` and introduce `concurrence`. ([#414], [#418], [#419])
- Introduce some metric functions. ([#414], [#420])
  - `hilbert_dist`
  - `hellinger_dist`
  - `bures_dist`
  - `bures_angle`
- Align `steadystate` ODE solver to other methods and improve GPU support. ([#421])

## [v0.27.0]
Release date: 2025-02-14

- Rename `sparse_to_dense` as `to_dense` and `dense_to_sparse` as `to_sparse`. ([#392])
- Fix erroneous definition of the stochastic term in `smesolve`. ([#393])
- Change name of `MultiSiteOperator` to `multisite_operator`. ([#394])
- Fix `smesolve` for specifying initial state as density matrix. ([#395])
- Add more generic solver for `steadystate_floquet` to allow more linear solvers. ([#396])
- Fix time evolution output when using `saveat` keyword argument. ([#398])
- Align some attributes of `mcsolve`, `ssesolve` and `smesolve` results with `QuTiP`. ([#402])
- Improve ensemble generation of `ssesolve` and change parameters handling on stochastic processes. ([#403])
- Set default trajectories to 500 and rename the keyword argument `ensemble_method` to `ensemblealg`. ([#405])
- Introduce measurement on `ssesolve` and `smesolve`. ([#404])

## [v0.26.0]
Release date: 2025-02-09

- Fix CUDA `sparse_to_dense`. ([#386])
- Improve pseudo inverse spectrum solver. ([#388])
- Add `smesolve` function for stochastic master equation. ([#389])

## [v0.25.2]
Release date: 2025-02-02

- Move code quality dependencies to separate environment. ([#380])
- Add additional normalization of the state during time evolution of `ssesolve`. This improves the numerical stability of the solver. ([#383])

## [v0.25.1]
Release date: 2025-01-29

- Fix Dynamical Fock Dimension states saving due to wrong saving of dimensions. ([#375])
- Support a list of observables for `expect`. ([#374], [#376])
- Add checks for `tlist` in time evolution solvers. The checks are to ensure that `tlist` is not empty, the elements are in increasing order, and the elements are unique. ([#378])

## [v0.25.0]
Release date: 2025-01-20

- Change the structure of block diagonalization functions, using `BlockDiagonalForm` struct and changing the function name from `bdf` to `block_diagonal_form`. ([#349])
- Add **GPUArrays** compatibility for `ptrace` function, by using **KernelAbstractions.jl**. ([#350])
- Introduce `Space`, `Dimensions`, `GeneralDimensions` structures to support wider definitions and operations of `Qobj/QobjEvo`, and potential functionalities in the future. ([#271], [#353], [#360])
- Improve lazy tensor warning for `SciMLOperators`. ([#370])
- Change order of `AbstractQuantumObject` data type. For example, from `QuantumObject{DataType,ObjType,DimsType}` to `QuantumObject{ObjType,DimsType,DataType}`. ([#371])

## [v0.24.0]
Release date: 2024-12-13

- Improve the construction of `QobjEvo`. ([#338], [#339])
- Support `Base.zero` and `Base.one` for `AbstractQuantumObject`. ([#342], [#346])
- Introduce visualization and function `plot_wigner` for easy plotting of Wigner functions. ([#86], [#292], [#347])

## [v0.23.1]
Release date: 2024-12-06

- Update `[compat]` to fix the incompatibility between `QuantumToolbox v0.22.0+` and `DiffEqCallbacks < v4.2.1`. ([#335])

## [v0.23.0]
Release date: 2024-12-04

- Change `SingleSiteOperator` with the more general `MultiSiteOperator`. ([#324])
- Make `spectrum` and `correlation` functions align with `Python QuTiP`, introduce spectrum solver `PseudoInverse`, remove spectrum solver `FFTCorrelation`, and introduce `spectrum_correlation_fft`. ([#330])

## [v0.22.0]
Release date: 2024-11-20

- Change the parameters structure of `sesolve`, `mesolve` and `mcsolve` functions to possibly support automatic differentiation. ([#311])
- Fix type instability and reduce extra memory allocation in `liouvillian`. ([#315], [#318])

## [v0.21.5]
Release date: 2024-11-15

- This is a demonstration of how to bump version number and also modify `CHANGELOG.md` before new release. ([#309])

## [v0.21.4]
Release date: 2024-11-13

- This is just a demonstration about [`Changelog.jl`](https://github.com/JuliaDocs/Changelog.jl). ([#139], [#306])


<!-- Links generated by Changelog.jl -->

[v0.21.4]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.4
[v0.21.5]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.21.5
[v0.22.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.22.0
[v0.23.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.0
[v0.23.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.23.1
[v0.24.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.24.0
[v0.25.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.0
[v0.25.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.1
[v0.25.2]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.25.2
[v0.26.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.26.0
[v0.27.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.27.0
[v0.28.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.28.0
[v0.29.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.29.0
[v0.29.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.29.1
[v0.30.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.30.0
[v0.30.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.30.1
[v0.31.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.31.0
[v0.31.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.31.1
[v0.32.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.32.0
[v0.32.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.32.1
[v0.33.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.33.0
[v0.34.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.34.0
[v0.34.1]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.34.1
[v0.35.0]: https://github.com/qutip/QuantumToolbox.jl/releases/tag/v0.35.0
[#86]: https://github.com/qutip/QuantumToolbox.jl/issues/86
[#139]: https://github.com/qutip/QuantumToolbox.jl/issues/139
[#271]: https://github.com/qutip/QuantumToolbox.jl/issues/271
[#292]: https://github.com/qutip/QuantumToolbox.jl/issues/292
[#306]: https://github.com/qutip/QuantumToolbox.jl/issues/306
[#309]: https://github.com/qutip/QuantumToolbox.jl/issues/309
[#311]: https://github.com/qutip/QuantumToolbox.jl/issues/311
[#315]: https://github.com/qutip/QuantumToolbox.jl/issues/315
[#318]: https://github.com/qutip/QuantumToolbox.jl/issues/318
[#324]: https://github.com/qutip/QuantumToolbox.jl/issues/324
[#330]: https://github.com/qutip/QuantumToolbox.jl/issues/330
[#335]: https://github.com/qutip/QuantumToolbox.jl/issues/335
[#338]: https://github.com/qutip/QuantumToolbox.jl/issues/338
[#339]: https://github.com/qutip/QuantumToolbox.jl/issues/339
[#342]: https://github.com/qutip/QuantumToolbox.jl/issues/342
[#346]: https://github.com/qutip/QuantumToolbox.jl/issues/346
[#347]: https://github.com/qutip/QuantumToolbox.jl/issues/347
[#349]: https://github.com/qutip/QuantumToolbox.jl/issues/349
[#350]: https://github.com/qutip/QuantumToolbox.jl/issues/350
[#353]: https://github.com/qutip/QuantumToolbox.jl/issues/353
[#360]: https://github.com/qutip/QuantumToolbox.jl/issues/360
[#370]: https://github.com/qutip/QuantumToolbox.jl/issues/370
[#371]: https://github.com/qutip/QuantumToolbox.jl/issues/371
[#374]: https://github.com/qutip/QuantumToolbox.jl/issues/374
[#375]: https://github.com/qutip/QuantumToolbox.jl/issues/375
[#376]: https://github.com/qutip/QuantumToolbox.jl/issues/376
[#378]: https://github.com/qutip/QuantumToolbox.jl/issues/378
[#380]: https://github.com/qutip/QuantumToolbox.jl/issues/380
[#383]: https://github.com/qutip/QuantumToolbox.jl/issues/383
[#386]: https://github.com/qutip/QuantumToolbox.jl/issues/386
[#388]: https://github.com/qutip/QuantumToolbox.jl/issues/388
[#389]: https://github.com/qutip/QuantumToolbox.jl/issues/389
[#392]: https://github.com/qutip/QuantumToolbox.jl/issues/392
[#393]: https://github.com/qutip/QuantumToolbox.jl/issues/393
[#394]: https://github.com/qutip/QuantumToolbox.jl/issues/394
[#395]: https://github.com/qutip/QuantumToolbox.jl/issues/395
[#396]: https://github.com/qutip/QuantumToolbox.jl/issues/396
[#398]: https://github.com/qutip/QuantumToolbox.jl/issues/398
[#402]: https://github.com/qutip/QuantumToolbox.jl/issues/402
[#403]: https://github.com/qutip/QuantumToolbox.jl/issues/403
[#404]: https://github.com/qutip/QuantumToolbox.jl/issues/404
[#405]: https://github.com/qutip/QuantumToolbox.jl/issues/405
[#408]: https://github.com/qutip/QuantumToolbox.jl/issues/408
[#410]: https://github.com/qutip/QuantumToolbox.jl/issues/410
[#411]: https://github.com/qutip/QuantumToolbox.jl/issues/411
[#413]: https://github.com/qutip/QuantumToolbox.jl/issues/413
[#414]: https://github.com/qutip/QuantumToolbox.jl/issues/414
[#416]: https://github.com/qutip/QuantumToolbox.jl/issues/416
[#418]: https://github.com/qutip/QuantumToolbox.jl/issues/418
[#419]: https://github.com/qutip/QuantumToolbox.jl/issues/419
[#420]: https://github.com/qutip/QuantumToolbox.jl/issues/420
[#421]: https://github.com/qutip/QuantumToolbox.jl/issues/421
[#423]: https://github.com/qutip/QuantumToolbox.jl/issues/423
[#428]: https://github.com/qutip/QuantumToolbox.jl/issues/428
[#430]: https://github.com/qutip/QuantumToolbox.jl/issues/430
[#436]: https://github.com/qutip/QuantumToolbox.jl/issues/436
[#437]: https://github.com/qutip/QuantumToolbox.jl/issues/437
[#438]: https://github.com/qutip/QuantumToolbox.jl/issues/438
[#440]: https://github.com/qutip/QuantumToolbox.jl/issues/440
[#443]: https://github.com/qutip/QuantumToolbox.jl/issues/443
[#448]: https://github.com/qutip/QuantumToolbox.jl/issues/448
[#450]: https://github.com/qutip/QuantumToolbox.jl/issues/450
[#453]: https://github.com/qutip/QuantumToolbox.jl/issues/453
[#455]: https://github.com/qutip/QuantumToolbox.jl/issues/455
[#456]: https://github.com/qutip/QuantumToolbox.jl/issues/456
[#460]: https://github.com/qutip/QuantumToolbox.jl/issues/460
[#470]: https://github.com/qutip/QuantumToolbox.jl/issues/470
[#472]: https://github.com/qutip/QuantumToolbox.jl/issues/472
[#473]: https://github.com/qutip/QuantumToolbox.jl/issues/473
[#476]: https://github.com/qutip/QuantumToolbox.jl/issues/476
[#480]: https://github.com/qutip/QuantumToolbox.jl/issues/480
[#485]: https://github.com/qutip/QuantumToolbox.jl/issues/485
[#486]: https://github.com/qutip/QuantumToolbox.jl/issues/486
[#487]: https://github.com/qutip/QuantumToolbox.jl/issues/487
[#489]: https://github.com/qutip/QuantumToolbox.jl/issues/489
[#494]: https://github.com/qutip/QuantumToolbox.jl/issues/494
[#500]: https://github.com/qutip/QuantumToolbox.jl/issues/500
[#501]: https://github.com/qutip/QuantumToolbox.jl/issues/501
[#504]: https://github.com/qutip/QuantumToolbox.jl/issues/504
[#506]: https://github.com/qutip/QuantumToolbox.jl/issues/506
[#507]: https://github.com/qutip/QuantumToolbox.jl/issues/507
[#509]: https://github.com/qutip/QuantumToolbox.jl/issues/509
[#512]: https://github.com/qutip/QuantumToolbox.jl/issues/512
[#513]: https://github.com/qutip/QuantumToolbox.jl/issues/513
[#515]: https://github.com/qutip/QuantumToolbox.jl/issues/515
[#517]: https://github.com/qutip/QuantumToolbox.jl/issues/517
[#520]: https://github.com/qutip/QuantumToolbox.jl/issues/520
[#531]: https://github.com/qutip/QuantumToolbox.jl/issues/531
[#536]: https://github.com/qutip/QuantumToolbox.jl/issues/536
[#537]: https://github.com/qutip/QuantumToolbox.jl/issues/537
[#539]: https://github.com/qutip/QuantumToolbox.jl/issues/539
