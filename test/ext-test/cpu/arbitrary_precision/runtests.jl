using QuantumToolbox
using ParallelTestRunner

QuantumToolbox.about()

testsuite = find_tests(dirname(@__FILE__))
delete!(testsuite, "setup")

runtests(QuantumToolbox, ARGS; testsuite)
