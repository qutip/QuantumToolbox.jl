using QuantumToolbox
using ParallelTestRunner

QuantumToolbox.about()

testsuite = find_tests(dirname(@__FILE__))
delete!(testsuite, "setup")

include(joinpath(@__DIR__, "..", "..", "..", "utils", "generate_test_worker.jl"))

runtests(QuantumToolbox, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
