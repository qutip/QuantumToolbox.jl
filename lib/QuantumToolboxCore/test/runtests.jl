using QuantumToolboxCore
using ParallelTestRunner

# make sure only QuantumToolboxCore is loaded for this test
(QuantumToolboxCore.QT_LIBRARIES != Module[QuantumToolboxCore]) &&
    error("This test should be individually testing the QuantumToolboxCore library. However, other libraries are loaded:\n$(QuantumToolboxCore.QT_LIBRARIES)")

QuantumToolboxCore.about()

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolboxCore, ARGS; testsuite)
