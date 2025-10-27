using Test
using TestItemRunner
using Pkg

const GROUP_LIST = String["All", "Core", "Code-Quality", "AutoDiff_Ext", "Makie_Ext", "CUDA_Ext", "Arbitrary-Precision"]

const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nThe allowed groups are: $GROUP_LIST\n"))

# Core tests
if (GROUP == "All") || (GROUP == "Core")
    import QuantumToolbox

    QuantumToolbox.about()

    println("\nStart running Core tests...\n")
    @run_package_tests verbose=true
end

########################################################################
# Use traditional Test.jl instead of TestItemRunner.jl for other tests #
########################################################################

const testdir = dirname(@__FILE__)

if (GROUP == "All") || (GROUP == "Code-Quality")
    Pkg.activate("core-test/code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    include(joinpath(testdir, "core-test", "code-quality", "code_quality.jl"))
end

if (GROUP == "AutoDiff_Ext")
    Pkg.activate("ext-test/cpu/autodiff")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    using QuantumToolbox
    using ForwardDiff
    using Zygote
    using Enzyme
    using SciMLSensitivity

    QuantumToolbox.about()

    include(joinpath(testdir, "ext-test", "cpu", "autodiff", "autodiff.jl"))
end

if (GROUP == "Makie_Ext")
    Pkg.activate("ext-test/cpu/makie")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    using QuantumToolbox
    QuantumToolbox.about()

    # CarioMakie is imported in the following script
    include(joinpath(testdir, "ext-test", "cpu", "makie", "makie_ext.jl"))
end

if (GROUP == "CUDA_Ext")
    Pkg.activate("ext-test/gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    using QuantumToolbox
    import LinearAlgebra: Diagonal
    import StaticArraysCore: SVector
    using CUDA
    using CUDA.CUSPARSE
    using CUDSS
    using LinearSolve

    QuantumToolbox.about()
    CUDA.versioninfo()

    include(joinpath(testdir, "ext-test", "gpu", "cuda_ext.jl"))
end

if (GROUP == "Arbitrary-Precision")
    Pkg.activate("ext-test/cpu/arbitrary_precision")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()

    using QuantumToolbox
    using SparseArrays
    using GenericSchur

    QuantumToolbox.about()

    include(joinpath(testdir, "ext-test", "cpu", "arbitrary_precision", "arbitrary_precision.jl"))
end
