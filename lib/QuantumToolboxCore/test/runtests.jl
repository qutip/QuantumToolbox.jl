using QuantumToolboxCore
using ParallelTestRunner

QuantumToolboxCore.about()

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolboxCore, ARGS; testsuite)
