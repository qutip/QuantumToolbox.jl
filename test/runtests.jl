using Test
using TestItemRunner
using Pkg

const testdir = dirname(@__FILE__)

# Define the paths to the library
const LIBRARY_NAME_AND_PATH = Dict(
    "Utils" => ("QuantumToolboxUtils", joinpath(testdir, "..", "lib", "QuantumToolboxUtils")),
)
const LIBRARY_LIST = collect(keys(LIBRARY_NAME_AND_PATH))

# Define the paths to the extension tests
const EXTENSION_PATH = Dict(
    "AutoDiff-Ext" => joinpath(testdir, "ext-test", "cpu", "autodiff"),
    "Makie-Ext" => joinpath(testdir, "ext-test", "cpu", "makie"),
    "CUDA-Ext" => joinpath(testdir, "ext-test", "gpu"),
    "Arbitrary-Precision" => joinpath(testdir, "ext-test", "cpu", "arbitrary_precision"),
)
const EXTENSION_LIST = collect(keys(EXTENSION_PATH))

# Handle the GROUP environment variable to determine which tests to run
const GROUP = get(ENV, "GROUP", "All")
const GROUP_LIST = String[
    "All",
    "Main",
    "Code-Quality",
    LIBRARY_LIST...,
    EXTENSION_LIST...,
]
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nThe allowed groups are: $GROUP_LIST\n"))

# function to set up the environment for subtests
function setup_subtest_env(path::String)
    Pkg.activate(path)
    if VERSION < v"1.11"
        for lib in LIBRARY_LIST
            _, lib_path = LIBRARY_NAME_AND_PATH[lib]
            Pkg.develop(path = lib_path)
        end
    end
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.update()
    return nothing
end

#######################################
# Main tests (use TestItemRunner.jl) #
#######################################
if (GROUP == "All") || (GROUP == "Main")
    import QuantumToolbox

    QuantumToolbox.about()

    println("\nStart running Main tests...")

    # tests in lib folder for each library
    # PATH: lib/LIBRARY_NAME/test/
    for lib in LIBRARY_LIST
        lib_name, _ = LIBRARY_NAME_AND_PATH[lib]
        println("\n[$lib_name]")
        @run_package_tests filter = ti -> occursin(joinpath("lib", lib_name, "test"), ti.filename) verbose = true
    end

    # main package tests (all tests except those in the lib folder)
    println("\n[QuantumToolbox]")
    @run_package_tests filter = ti -> !occursin("lib", ti.filename) verbose = true

    println("\n===============> Main tests completed <===============\n")
end

# only run tests for a specific library
if GROUP ∈ LIBRARY_LIST
    import QuantumToolbox

    QuantumToolbox.about()

    lib_name, _ = LIBRARY_NAME_AND_PATH[GROUP]
    println("\n[$lib_name]")
    @run_package_tests filter = ti -> occursin(joinpath("lib", lib_name, "test"), ti.filename) verbose = true
end

########################################################################
# Use traditional Test.jl instead of TestItemRunner.jl for other tests #
########################################################################

# Code Quality tests
if (GROUP == "All") || (GROUP == "Code-Quality")
    path = joinpath(testdir, "code-quality")
    setup_subtest_env(path)

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    # run code quality tests for all libraries and the main package
    for lib in LIBRARY_LIST
        _, lib_path = LIBRARY_NAME_AND_PATH[lib]
        include(joinpath(lib_path, "test", "code_quality.jl"))
    end
    include(joinpath(path, "code_quality.jl"))
end

# Extension tests
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATH[GROUP]
    setup_subtest_env(path)

    if GROUP == "AutoDiff-Ext"
        println(Pkg.status())
        include(joinpath(path, "autodiff.jl"))
    elseif GROUP == "Makie-Ext"
        include(joinpath(path, "makie_ext.jl"))
    elseif GROUP == "CUDA-Ext"
        include(joinpath(path, "cuda_ext.jl"))
    elseif GROUP == "Arbitrary-Precision"
        include(joinpath(path, "arbitrary_precision.jl"))
    end
end
