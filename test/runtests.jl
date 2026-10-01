using ParallelTestRunner
using Pkg

# include GROUP_LIST and test paths
include("utils/group_list.jl")

# Handle the GROUP environment variable to determine which tests to run
const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nAvailable test GROUP are:\n$SHOW_GROUP_LIST\n"))

# function to set up the environment for subtests
function setup_subtest_env(path::String)
    Pkg.activate(path)

    # `cd` into `path` first; otherwise `Pkg.develop` writes an absolute (machine-specific) path to Project.toml
    cd(path) do
        Pkg.develop(PackageSpec(path = relpath(dirname(@__DIR__), path))) # must `develop` otherwise the code coverage for libraries (like CUDA) will not be collected
    end
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

    # Main package tests
    main_tests = find_tests(joinpath(testdir, "main-test"))
    delete!(main_tests, "basic_solvers/setup") # remove the setup.jl, since it is not a test but a setup file

    for (name, include_expr) in main_tests
        testsuite["Main/$name"] = include_expr
    end

    # tests in lib folder for each library
    # PATH: lib/***/test/
    for lib in LIBRARY_LIST
        path = LIBRARY_PATH[lib]
        lib_tests = find_tests(path)
        delete!(lib_tests, "setup") # remove the setup.jl, since it is not a test but a setup file

        for (name, include_expr) in lib_tests
            testsuite["$lib/$name"] = include_expr
        end
    end

    import QuantumToolbox
    QuantumToolbox.about()
    println("[Tests for GROUP = $GROUP]")

    include("utils/generate_test_worker.jl")
    runtests(QuantumToolbox, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
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

    println("[Tests for GROUP = $GROUP]")
    include(joinpath(path, "code_quality.jl"))
end

##############################
# (individual) Library tests #
##############################
if GROUP ∈ LIBRARY_LIST
    lib_path = LIBRARY_PATH[GROUP]
    setup_subtest_env(lib_path)

    println("[Tests for GROUP = $GROUP]")
    include(joinpath(lib_path, "runtests.jl"))
end

###################
# Extension tests #
###################
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATH[GROUP]
    setup_subtest_env(path)

    println("[Tests for GROUP = $GROUP]")
    include(joinpath(path, "runtests.jl"))
end
