using ParallelTestRunner
using Pkg

const testdir = dirname(@__FILE__)

# Define the paths to the library-tests
const LIBRARY_PATH = Dict(
    "Core" => joinpath(testdir, "..", "lib", "QuantumToolboxCore", "test"),
    "Visual" => joinpath(testdir, "..", "lib", "QuantumToolboxVisual", "test"),
)
const LIBRARY_LIST = collect(keys(LIBRARY_PATH))

# Define the paths to the extension tests
const EXTENSION_PATH = Dict(
    "AutoDiff-Ext" => joinpath(testdir, "ext-test", "cpu", "autodiff"),
    "Makie-Ext" => joinpath(testdir, "ext-test", "cpu", "makie"),
    "CUDA-Ext" => joinpath(testdir, "ext-test", "cuda"),
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
    Pkg.update()
    return nothing
end

######################################
# Main package and all library tests #
######################################
if (GROUP == "All") || (GROUP == "Main")
    # build up the set of tests to run for this GROUP, merging the main package's
    # tests with all library-tests (namespaced by library name, e.g., "Core/quantum_objects")
    testsuite = Dict{String, Expr}()

    if (GROUP == "All") || (GROUP == "Main")
        main_tests = find_tests(joinpath(testdir, "main-test"))

        for (name, include_expr) in main_tests
            testsuite["Main/$name"] = include_expr
        end
    end

    # tests in lib folder for each library
    # PATH: lib/***/test/
    for lib in LIBRARY_LIST
        path = LIBRARY_PATH[lib]
        lib_tests = find_tests(path)

        for (name, include_expr) in lib_tests
            testsuite["$lib/$name"] = include_expr
        end
    end

    import QuantumToolbox
    QuantumToolbox.about()
    runtests(QuantumToolbox, ARGS; testsuite)
end

######################
# Code Quality tests #
######################
if (GROUP == "All") || (GROUP == "Code-Quality")
    path = joinpath(testdir, "code-quality")
    setup_subtest_env(path)

    using QuantumToolbox
    using Aqua, JET

    (GROUP == "Code-Quality") && QuantumToolbox.about() # print version info. for code quality CI in GitHub

    include(joinpath(path, "code_quality.jl"))
end

##############################
# (individual) Library tests #
##############################
if GROUP ∈ LIBRARY_LIST
    lib_path = LIBRARY_PATH[GROUP]
    setup_subtest_env(lib_path)

    include(joinpath(lib_path, "runtests.jl"))
end

###################
# Extension tests #
###################
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATH[GROUP]
    setup_subtest_env(path)

    include(joinpath(path, "runtests.jl"))
end
