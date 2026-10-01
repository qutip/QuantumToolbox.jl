using QuantumToolboxCore
using ParallelTestRunner

# make sure only QuantumToolboxCore is loaded for this test
(QuantumToolboxCore.QT_LIBRARIES != Module[QuantumToolboxCore]) &&
    error("This test should be individually testing the QuantumToolboxCore library. However, other libraries are loaded:\n$(QuantumToolboxCore.QT_LIBRARIES)")

QuantumToolboxCore.about()

testsuite = find_tests(dirname(@__FILE__))

include(joinpath(@__DIR__, "..", "..", "..", "test", "utils", "generate_test_worker.jl"))

runtests(QuantumToolboxCore, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
