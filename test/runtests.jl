using Test
using Pkg

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

# Put core tests in alphabetical order
core_tests = [
    "block_diagonal_form.jl",
    "correlations_and_spectrum.jl",
    "dynamical_fock_dimension_mesolve.jl", # TODO: fix tests
    "dynamical-shifted-fock.jl",
    "eigenvalues_and_operators.jl",
    "entanglement.jl",
    "generalized_master_equation.jl",
    "low_rank_dynamics.jl",
    "negativity_and_partial_transpose.jl",
    "progress_bar.jl",
    "quantum_objects.jl",
    "quantum_objects_evo.jl",
    "states_and_operators.jl",
    "steady_state.jl",
    "time_evolution.jl",
    "utilities.jl",
    "wigner.jl",
]

if (GROUP == "All") || (GROUP == "Code-Quality")
    using QuantumToolbox
    using Aqua, JET

    include(joinpath(testdir, "core-test", "code_quality.jl"))
end

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

if (GROUP == "All") || (GROUP == "CairoMakie_Ext")
    using QuantumToolbox

    (GROUP == "CairoMakie_Ext") && QuantumToolbox.about()

    # CarioMakie is imported in the following script
    include(joinpath(testdir, "ext-test", "cairomakie_ext.jl"))
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
