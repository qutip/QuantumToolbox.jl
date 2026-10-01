using QuantumToolboxVisual
using ParallelTestRunner

# make sure only QuantumToolboxVisual is loaded for this test
(QuantumToolboxCore.QT_LIBRARIES != Module[QuantumToolboxVisual, QuantumToolboxCore]) &&
    error("This test should be testing the QuantumToolboxVisual library individually. However, extra libraries are loaded:\n$(QuantumToolboxCore.QT_LIBRARIES)")

QuantumToolboxVisual.about()

testsuite = find_tests(dirname(@__FILE__))

include(joinpath(@__DIR__, "..", "..", "..", "test", "utils", "generate_test_worker.jl"))

runtests(QuantumToolboxVisual, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
