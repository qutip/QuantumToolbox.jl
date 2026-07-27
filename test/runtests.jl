using Test
using TestItemRunner
using Pkg

const testdir = dirname(@__FILE__)

# Define the paths to the extension tests
const EXTENSION_PATHS = Dict(
    "AutoDiff_Ext" => joinpath(testdir, "ext-test", "cpu", "autodiff"),
    "Makie_Ext" => joinpath(testdir, "ext-test", "cpu", "makie"),
    "CUDA_Ext" => joinpath(testdir, "ext-test", "gpu"),
    "Arbitrary-Precision" => joinpath(testdir, "ext-test", "cpu", "arbitrary_precision"),
)
const EXTENSION_LIST = collect(keys(EXTENSION_PATHS))

# Handle the GROUP environment variable to determine which tests to run
const GROUP = get(ENV, "GROUP", "All")
const GROUP_LIST = String[
    "All",
    "Core",
    "Code-Quality",
    EXTENSION_LIST...,
]
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nThe allowed groups are: $GROUP_LIST\n"))

# function to set up the environment for subtests
function setup_subtest_env(path::String)
    Pkg.activate(path)
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()
    return nothing
end

######################################
# Core tests (use TestItemRunner.jl) #
######################################
if (GROUP == "All") || (GROUP == "Core")
    import QuantumToolbox

    QuantumToolbox.about()

    println("\nStart running Core tests...\n")
    @run_package_tests verbose = true
end

########################################################################
# Use traditional Test.jl instead of TestItemRunner.jl for other tests #
########################################################################

# Code Quality tests
if (GROUP == "All") || (GROUP == "Code-Quality")
    path = joinpath(testdir, "core-test", "code-quality")
    setup_subtest_env(path)

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    include(joinpath(path, "code_quality.jl"))
end

# Extension tests
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATHS[GROUP]
    setup_subtest_env(path)

    if GROUP == "AutoDiff_Ext"
        println(Pkg.status())
        include(joinpath(path, "autodiff.jl"))
    elseif GROUP == "Makie_Ext"
        include(joinpath(path, "makie_ext.jl"))
    elseif GROUP == "CUDA_Ext"
        include(joinpath(path, "cuda_ext.jl"))
    elseif GROUP == "Arbitrary-Precision"
        include(joinpath(path, "arbitrary_precision.jl"))
    end
end
