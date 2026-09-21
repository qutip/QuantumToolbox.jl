using QuantumToolboxVisual
using ParallelTestRunner

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolboxVisual, ARGS; testsuite)
