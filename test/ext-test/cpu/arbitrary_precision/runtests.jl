using QuantumToolbox
using ParallelTestRunner

QuantumToolbox.about()

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolbox, ARGS; testsuite)
