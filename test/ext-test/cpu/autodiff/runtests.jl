using QuantumToolbox
using ParallelTestRunner
import Pkg

println(Pkg.status())

QuantumToolbox.about()

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolbox, ARGS; testsuite)
