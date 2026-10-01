using QuantumToolboxVisual
using ParallelTestRunner

testsuite = find_tests(dirname(@__FILE__))

include(joinpath(@__DIR__, "..", "..", "..", "utils", "generate_test_worker.jl"))

runtests(QuantumToolboxVisual, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
