using Test
using ParallelTestRunner
using Pkg

const testdir = dirname(@__FILE__)

# Define the paths to the library
const LIBRARY_NAME_AND_PATH = Dict(
    "Core" => ("QuantumToolboxCore", joinpath(testdir, "..", "lib", "QuantumToolboxCore")),
    "Visual" => ("QuantumToolboxVisual", joinpath(testdir, "..", "lib", "QuantumToolboxVisual")),
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
    specs = PackageSpec[]
    if VERSION < v"1.11"
        for lib in LIBRARY_LIST
            _, lib_path = LIBRARY_NAME_AND_PATH[lib]
            push!(specs, PackageSpec(path = lib_path))
        end
    end
    push!(specs, PackageSpec(path = dirname(@__DIR__)))
    Pkg.develop(specs)
    Pkg.update()
    return nothing
end

##################################
# Main package and library tests #
##################################
if (GROUP == "All") || (GROUP == "Main") || (GROUP ∈ LIBRARY_LIST)
    # build up the set of tests to run for this GROUP, merging the main package's
    # tests with each requested library's tests (namespaced by library name, e.g., "Core/quantum_objects")
    testsuite = Dict{String, Expr}()

    if (GROUP == "All") || (GROUP == "Main")
        main_tests = find_tests(joinpath(testdir, "main-test"))

        for (name, include_expr) in main_tests
            testsuite["Main/$name"] = include_expr
        end
    end

    # tests in lib folder for each library
    # PATH: lib/LIBRARY_NAME/test/
    for lib in LIBRARY_LIST
        (GROUP == "All") || (GROUP == "Main") || (GROUP == lib) || continue

        lib_name, lib_path = LIBRARY_NAME_AND_PATH[lib]

        lib_tests = find_tests(joinpath(lib_path, "test"))

        for (name, include_expr) in lib_tests
            testsuite["$lib/$name"] = include_expr
        end
    end

    # only import the package that `runtests` needs as its module argument (used
    # solely to name the historical-duration cache file, never `using`-ed into any
    # worker) -- a `Core`-only or `Visual`-only run must not load the main
    # QuantumToolbox package at all, in the coordinator process or in any worker
    if GROUP == "Core"
        import QuantumToolboxCore
        QuantumToolboxCore.about()
        runtests(QuantumToolboxCore, ARGS; testsuite)
    elseif GROUP == "Visual"
        import QuantumToolboxVisual
        QuantumToolboxVisual.about()
        runtests(QuantumToolboxVisual, ARGS; testsuite)
    else # "All" or "Main"
        import QuantumToolbox
        QuantumToolbox.about()
        runtests(QuantumToolbox, ARGS; testsuite)
    end
end

############################################################
# Use traditional Test.jl instead of ParallelTestRunner.jl #
############################################################

# Code Quality tests
if (GROUP == "All") || (GROUP == "Code-Quality")
    path = joinpath(testdir, "code-quality")
    setup_subtest_env(path)

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    include(joinpath(path, "code_quality.jl"))
end

###################
# Extension tests #
###################
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATH[GROUP]
    setup_subtest_env(path)

    (GROUP == "AutoDiff-Ext") && println(Pkg.status())

    include(joinpath(path, "runtests.jl"))
end
