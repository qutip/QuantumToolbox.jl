using Test
using Pkg

# Importing only the necessary functions to keep track the re-export of the functions
import LinearAlgebra: Diagonal, I, mul!, triu, tril, triu!, tril!
import SparseArrays: sparse, sprand, spzeros, spdiagm, nnz, SparseVector, SparseMatrixCSC, AbstractSparseMatrix
import StaticArraysCore: SVector

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

# Put core tests in alphabetical order
core_tests = [
    # "block_diagonal_form.jl",
    # "correlations_and_spectrum.jl",
    # "dynamical_fock_dimension_mesolve.jl",
    # "dynamical-shifted-fock.jl",
    # "eigenvalues_and_operators.jl",
    # "entropy_and_metric.jl",
    # "generalized_master_equation.jl",
    # "low_rank_dynamics.jl",
    # "negativity_and_partial_transpose.jl",
    # "progress_bar.jl",
    # "quantum_objects.jl",
    # "quantum_objects_evo.jl",
    # "states_and_operators.jl",
    # "steady_state.jl",
    "time_evolution.jl",
    # "utilities.jl",
    # "wigner.jl",
]

if (GROUP == "All") || (GROUP == "Core")
    using QuantumToolbox
    import QuantumToolbox: position, momentum
    import Random: MersenneTwister
    import SciMLOperators: MatrixOperator, NullOperator, IdentityOperator

    QuantumToolbox.about()

    for test in core_tests
        include(joinpath(testdir, "core-test", test))
    end
end

if (GROUP == "All") || (GROUP == "Code-Quality")
    Pkg.activate("core-test/code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using QuantumToolbox
    using Aqua, JET

    include(joinpath(testdir, "core-test", "code-quality", "code_quality.jl"))
end

if (GROUP == "AutoDiff_Ext")
    Pkg.activate("ext-test/cpu/autodiff")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using QuantumToolbox
    using Zygote
    using Enzyme
    using SciMLSensitivity

    include(joinpath(testdir, "ext-test", "cpu", "autodiff", "zygote.jl"))
end

if (GROUP == "Makie_Ext")
    Pkg.activate("ext-test/cpu/makie")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using QuantumToolbox
    QuantumToolbox.about()

    # CarioMakie is imported in the following script
    include(joinpath(testdir, "ext-test", "cpu", "makie", "makie_ext.jl"))
end

if (GROUP == "CUDA_Ext")
    Pkg.activate("ext-test/gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using QuantumToolbox
    using CUDA
    using CUDA.CUSPARSE
    # CUDA.allowscalar(false) # This is already set in the extension script

    QuantumToolbox.about()
    CUDA.versioninfo()

    include(joinpath(testdir, "ext-test", "gpu", "cuda_ext.jl"))
end
