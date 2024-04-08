using Test
using Pkg
using QuantumToolbox

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
    "steady_state.jl",
    "time_evolution_and_partial_trace.jl",
    "wigner.jl",
]

if (GROUP == "All") || (GROUP == "Core")
    for test in core_tests
        include(joinpath(testdir, test))
    end
end

if (GROUP == "All") || (GROUP == "Code-Quality")
    Pkg.activate("code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include(joinpath(testdir, "code-quality/aqua.jl"))
    include(joinpath(testdir, "code-quality/jet.jl"))
end