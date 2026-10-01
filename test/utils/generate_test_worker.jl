import ParallelTestRunner: addworker

# marker comment that a test file can include (anywhere in the file) to use customized workers
const MARKER_MULTITHREAD = "#!PTR_MULTITHREAD"

function generate_test_worker(testsuite::Dict{String, Expr})
    # find out all test files that contain #!PTR_MULTITHREAD
    multithread_tests = Set{String}()
    for (name, ex) in testsuite
        path = ex.args[2]
        if any(line -> strip(line) == MARKER_MULTITHREAD, eachline(path))
            push!(multithread_tests, name)
        end
    end

    # return a test_worker function
    return function (name::String)
        if name in multithread_tests
            return addworker(; exeflags = ["--threads=4"])
        else
            return nothing # uses default worker
        end
    end
end
