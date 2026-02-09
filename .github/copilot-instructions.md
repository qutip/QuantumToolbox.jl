# QuantumToolbox.jl - Copilot Instructions

## Project Overview

QuantumToolbox.jl is a Julia package for quantum physics simulations, mirroring Python's QuTiP. It provides quantum state/operator manipulation, solvers for time evolution (Schrödinger, Lindblad, Monte Carlo, stochastic), GPU support, distributed computing, and automatic differentiation via SciML ecosystem.

## Architecture Essentials

### Core Type System: Quantum Objects

The project centers on `AbstractQuantumObject` with two main concrete types:
- **`QuantumObject`**: Time-independent quantum objects (states, operators, superoperators)
  - Located: [src/qobj/quantum_object.jl](../src/qobj/quantum_object.jl)
  - Core attributes: `data` (AbstractArray), `type` (QuantumObjectType), `dimensions` (Dimensions)
  - Type tags: `Ket()`, `Bra()`, `Operator()`, `SuperOperator()`, `OperatorKet()`, `OperatorBra()`
  
- **`QuantumObjectEvolution`** (QobjEvo): Time-dependent quantum objects
  - Located: [src/qobj/quantum_object_evo.jl](../src/qobj/quantum_object_evo.jl)
  - Wraps SciMLOperators for efficient coefficient updates

**Key principle**: Type tags are ALWAYS required. Never create QuantumObjects without specifying type and dimensions.

### Solver Architecture

Time evolution solvers follow SciML pattern: `*Problem` → `solve()` → `*Solution`
- **Master Equation (mesolve)**: Lindbladian + collapse operators
- **Schrödinger (sesolve)**: Pure state evolution  
- **Monte Carlo (mcsolve)**: Stochastic quantum trajectories
- **Stochastic (smesolve, ssesolve)**: Brownian/Wiener noise
- **Bloch-Redfield (brmesolve)**: Spectral environment coupling

All located in [src/time_evolution/](../src/time_evolution/). Each implements same pattern:
```julia
prob = mesolveProblem(H, ψ0, tlist, c_ops; ...)  # Create problem
sol = mesolve(H, ψ0, tlist, c_ops; ...)          # Solve directly
```

### Module Structure

- **[src/](../src/)**: Core quantum computing logic
  - `qobj/`: Type definitions and object operations
  - `time_evolution/`: Solver implementations
  - `utilities.jl`: Helper functions (tensor products, partial trace, etc.)
  - `steadystate.jl`, `entropy.jl`, `correlations.jl`: Analysis tools
  
- **[ext/](../ext/)**: Optional extensions (loaded via weakdeps)
  - `QuantumToolboxMakieExt.jl`: Visualization (Bloch sphere, Wigner)
  - `QuantumToolboxCUDAExt.jl`: GPU acceleration
  - `QuantumToolboxChainRulesCoreExt.jl`: Automatic differentiation

- **[test/](../test/)**: Two-tier test structure
  - `core-test/`: Core functionality via `TestItemRunner.jl`
  - `ext-test/`: Extensions (cpu/, gpu/ subdirs)

## Code Patterns & Conventions

### Sparse Matrix Usage
- Default data storage: `SparseMatrixCSC{ComplexF64}` or `Vector{ComplexF64}`
- Use `SparseArrays.jl` functions: `nnz()`, `nonzeros()`, `droptol!()`, `dropzeros!()`
- Auto-tidyup controlled by `QuantumToolbox.settings.auto_tidyup` (default: true, tolerance: 1e-14)

### Dimension Handling
Never assume 1D systems. All code must handle arbitrary composite Hilbert spaces:
- `dimensions::Dimensions` stores as `StaticArraysCore.SVector`
- Access via `qobj.dims` (returns dims Tuple) or `qobj.dimensions` (full Dimensions object)
- Use `Dimensions()` constructor to create from integer/tuple specifications

### Error Checking
Use `@assert` for internal checks; throw descriptive `ArgumentError` for invalid user inputs. Example:
```julia
_check_QuantumObject(type, dimensions, rows, cols)  # In quantum_object.jl
```

### Progress Bars
Implemented via `ProgressMeter.jl` with settings at [src/settings.jl](../src/settings.jl):
```julia
QuantumToolbox.settings.ProgressMeterKWARGS = (showspeed=true, printed=true)
```

## Developer Workflows

### Running Tests
```bash
make test                    # Run all core tests
GROUP=Core make test         # Run core tests only
GROUP=CUDA_Ext make test     # Run CUDA extension tests only
```

Test groups (set `GROUP` env var): "All", "Core", "Code-Quality", "AutoDiff_Ext", "Makie_Ext", "CUDA_Ext", "Arbitrary-Precision"

### Code Quality & Formatting
```bash
make format                  # Format code with Runic
make setup                   # Install dev dependencies (Runic, Changelog)
```

Uses Aqua.jl and JET.jl for static analysis (see [test/core-test/code-quality/](../test/core-test/code-quality/)).

### Documentation
```bash
make docs                    # Build docs (updates dependencies, runs make.jl)
make vitepress               # Start live preview (http://localhost:5173)
```

Documentation uses Documenter.jl + DocumenterVitepress + CairoMakie for examples. Doctests enabled by default.

### Benchmarks
Located in [benchmarks/](../benchmarks/):
```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl
```

## Integration Points & Dependencies

### SciML Integration
- **SciMLBase**: Problem/solution types, ODE/SDE functions, ensemble algorithms
- **OrdinaryDiffEq**: Vern7, DP5 algorithms (imported directly)
- **LinearSolve**: Krylov methods (KrylovJL_MINRES, GMRES) for steady-state
- **SciMLOperators**: For QobjEvo coefficient caching and updates

### External Packages
- **FFTW.jl**: FFT for spectrum/correlation functions
- **Graphs.jl**: For spin lattice operations
- **DifferentialEquations.jl**: StochasticDiffEq, DiffEqCallbacks

### Extension Dependencies
- CUDA.jl, GPUArrays.jl, KernelAbstractions.jl → GPU arrays
- Makie.jl → 3D visualization
- ChainRulesCore.jl → AD support

## Key Files to Reference

- **Type system**: [src/qobj/quantum_object_base.jl](../src/qobj/quantum_object_base.jl) (QuantumObjectType hierarchy)
- **Main type**: [src/qobj/quantum_object.jl](../src/qobj/quantum_object.jl) (QuantumObject struct)
- **Time-dependent**: [src/qobj/quantum_object_evo.jl](../src/qobj/quantum_object_evo.jl) (QobjEvo wrapper)
- **Solver examples**: [src/time_evolution/mesolve.jl](../src/time_evolution/mesolve.jl) (pattern template)
- **Settings**: [src/settings.jl](../src/settings.jl) (global configuration)
- **Main module**: [src/QuantumToolbox.jl](../src/QuantumToolbox.jl) (exports, includes)

## Common Tasks

**Adding new solver**: Copy mesolve (or similars) pattern, create `*Problem`/`solve` functions in new file, add to [src/QuantumToolbox.jl](../src/QuantumToolbox.jl) include list.

**Adding GPU support**: Extend ext/QuantumToolboxCUDAExt.jl with appropriate CuArray dispatches.

**Debugging**: When debugging, create a new environment in a new folder with only QuantumToolbox.jl and dependencies to isolate issues. Remove them when debugging is complete.

**Type-stable code**: Always use concrete types. Test with JET.jl. Use `@code_warntype` for debugging.

**Documentation**: Add docstrings with `@doc raw"""..."""` blocks in function definitions; examples auto-tested.

**Auto updating instructions**: Keep this file updated with any architectural or workflow changes.
