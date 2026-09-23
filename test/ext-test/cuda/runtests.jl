using QuantumToolbox
using CUDA
using ParallelTestRunner

QuantumToolbox.about()
CUDA.versioninfo()

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolbox, ARGS; testsuite, init_worker_code)
