using Test
using TestItemRunner
using Pkg

const EXTENSIONS = String[
    "AutoDiff_Ext",
    "Makie_Ext",
    "CUDA_Ext",
    "Arbitrary-Precision"
]

const GROUP_LIST = String[
    "All",
    "Core",
    "Code-Quality",
    EXTENSIONS...,
]

const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nThe allowed groups are: $GROUP_LIST\n"))

function setup_subtest_env(path::String)
    Pkg.activate(path)
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()
    return nothing
end

# Core tests
if (GROUP == "All") || (GROUP == "Core")
    import QuantumToolbox

    QuantumToolbox.about()

    println("\nStart running Core tests...\n")
    @run_package_tests verbose = true
end

########################################################################
# Use traditional Test.jl instead of TestItemRunner.jl for other tests #
########################################################################

const testdir = dirname(@__FILE__)

if (GROUP == "All") || (GROUP == "Code-Quality")
    path = "core-test/code-quality"
    setup_subtest_env()

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    include(joinpath(testdir, path, "code_quality.jl"))
end

if (GROUP == "AutoDiff_Ext")
    path = "ext-test/cpu/autodiff"
    setup_subtest_env(path)
    println(Pkg.status())

    using QuantumToolbox
    using ForwardDiff
    using Enzyme
    using Mooncake
    using SciMLSensitivity
    using SciMLSensitivity: MooncakeVJP

    QuantumToolbox.about()

    include(joinpath(testdir, path, "autodiff.jl"))
end

if (GROUP == "Makie_Ext")
    path = "ext-test/cpu/makie"
    setup_subtest_env(path)

    using QuantumToolbox
    QuantumToolbox.about()

    # CarioMakie is imported in the following script
    include(joinpath(testdir, path, "makie_ext.jl"))
end

if (GROUP == "CUDA_Ext")
    path = "ext-test/gpu"
    setup_subtest_env(path)

    using QuantumToolbox
    import LinearAlgebra: Diagonal
    import SparseArrays: SparseMatrixCSC
    using CUDA
    using CUDA.CUSPARSE
    using CUDSS
    using LinearSolve

    QuantumToolbox.about()
    CUDA.versioninfo()

    include(joinpath(testdir, path, "cuda_ext.jl"))
end

if (GROUP == "Arbitrary-Precision")
    path = "ext-test/cpu/arbitrary_precision"
    setup_subtest_env(path)

    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    setprecision(128) # Instead of 256. This speeds up the tests.

    using QuantumToolbox
    using LinearAlgebra
    using SparseArrays
    using Sparspak
    using GenericSchur

    QuantumToolbox.about()

    include(joinpath(testdir, path, "arbitrary_precision.jl"))
end
