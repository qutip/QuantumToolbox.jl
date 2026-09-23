const testdir = dirname(@__FILE__)

# Define the paths to the library-tests
const LIBRARY_PATH = Dict(
    "Core" => joinpath(testdir, "..", "lib", "QuantumToolboxCore", "test"),
    "Visual" => joinpath(testdir, "..", "lib", "QuantumToolboxVisual", "test"),
)
const LIBRARY_LIST = collect(keys(LIBRARY_PATH))

# Define the paths to the extension tests
const EXTENSION_PATH = Dict(
    "AutoDiff" => joinpath(testdir, "ext-test", "cpu", "autodiff"),
    "Makie" => joinpath(testdir, "ext-test", "cpu", "makie"),
    "CUDA" => joinpath(testdir, "ext-test", "cuda"),
    "Arbitrary-Precision" => joinpath(testdir, "ext-test", "cpu", "arbitrary_precision"),
)
const EXTENSION_LIST = collect(keys(EXTENSION_PATH))

# the list of test groups
const GROUP_LIST = String[
    "All",
    "Main",
    "Code-Quality",
    LIBRARY_LIST...,
    EXTENSION_LIST...,
]

# prettier display of the list of test groups
const SHOW_GROUP_LIST = join("- " .* GROUP_LIST, "\n")

# if this file is executed directly, print the available test groups
if abspath(PROGRAM_FILE) == @__FILE__
    println("Available test GROUP:\n$(SHOW_GROUP_LIST)")
    nothing
end
