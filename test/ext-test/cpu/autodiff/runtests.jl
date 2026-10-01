using QuantumToolbox
using ParallelTestRunner
import Pkg

println(Pkg.status())

QuantumToolbox.about()

testsuite = find_tests(dirname(@__FILE__))

include(joinpath(@__DIR__, "..", "..", "..", "utils", "generate_test_worker.jl"))

runtests(QuantumToolbox, ARGS; testsuite, test_worker = generate_test_worker(testsuite))
