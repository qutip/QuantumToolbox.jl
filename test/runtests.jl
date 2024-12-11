using Test
using Pkg

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

# Put core tests in alphabetical order
core_tests = [
    "correlations_and_spectrum.jl",
    "dynamical_fock_dimension_mesolve.jl",
    "dynamical-shifted-fock.jl",
    "eigenvalues_and_operators.jl",
    "entanglement.jl",
    "generalized_master_equation.jl",
    "low_rank_dynamics.jl",
    "negativity_and_partial_transpose.jl",
    "permutation.jl",
    "progress_bar.jl",
    "quantum_objects.jl",
    "quantum_objects_evo.jl",
    "states_and_operators.jl",
    "steady_state.jl",
    "time_evolution.jl",
    "utilities.jl",
    "wigner.jl",
]

ext_tests = [joinpath("cairomakie", "cairomakie_ext.jl")]

if (GROUP == "All") || (GROUP == "Code-Quality")
    using QuantumToolbox
    using Aqua, JET

    include(joinpath(testdir, "core-test", "code_quality.jl"))
end

if (GROUP == "All") || (GROUP == "Core")
    using QuantumToolbox
    import QuantumToolbox: position, momentum
    import Random: MersenneTwister
    import SciMLOperators: MatrixOperator

    QuantumToolbox.about()

    for test in core_tests
        include(joinpath(testdir, "core-test", test))
    end

    for test in ext_tests
        include(joinpath(testdir, "ext-test", test))
    end
end

if (GROUP == "CUDA_Ext")# || (GROUP == "All")
    Pkg.activate("ext-test/gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using QuantumToolbox
    using CUDA
    using CUDA.CUSPARSE
    CUDA.allowscalar(false) # Avoid unexpected scalar indexing

    QuantumToolbox.about()
    CUDA.versioninfo()

    include(joinpath(testdir, "ext-test", "gpu", "cuda_ext.jl"))
end
