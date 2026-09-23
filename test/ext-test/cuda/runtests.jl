using QuantumToolbox
using CUDA
using ParallelTestRunner

QuantumToolbox.about()
CUDA.versioninfo()

# TEMPORARY DIAGNOSTIC: confirm whether --code-coverage propagates from this
# coordinator process to the Malt-spawned workers, and whether the paths
# involved actually agree. Remove once resolved.
function _coverage_diagnostic(label)
    opts = Base.JLOptions()
    tp = opts.tracked_path == C_NULL ? "<none>" : unsafe_string(opts.tracked_path)
    @info "$label coverage opts" code_coverage=opts.code_coverage tracked_path=tp pwd=pwd() this_file=@__FILE__
end
_coverage_diagnostic("coordinator")
const init_worker_code = quote
    opts = Base.JLOptions()
    tp = opts.tracked_path == C_NULL ? "<none>" : unsafe_string(opts.tracked_path)
    @info "worker coverage opts" code_coverage=opts.code_coverage tracked_path=tp pwd=pwd() this_file=@__FILE__
end

testsuite = find_tests(dirname(@__FILE__))

runtests(QuantumToolbox, ARGS; testsuite, init_worker_code)
