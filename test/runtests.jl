using TestItemRunner
using Pkg

using QuantumToolbox

const GROUP_LIST = String["All", "Core", "Code-Quality", "AutoDiff_Ext", "Makie_Ext", "CUDA_Ext"]

const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP"))

testfilter = ti -> begin
    if (GROUP == "All") || (GROUP == "Core")
        return :core in ti.tags
    end

    if GROUP == "AutoDiff_Ext"
        return :autodiff in ti.tags
    end

    if GROUP == "Makie_Ext"
        return :makie in ti.tags
    end

    if GROUP == "CUDA_Ext"
        return :cuda in ti.tags
    end

    if GROUP == "Code-Quality"
        return false
    end
end

QuantumToolbox.about()

println("\nStart running tests [for GROUP = $GROUP]...\n")

# TestItemRunner.jl
@run_package_tests filter=testfilter
(GROUP == "Code-Quality") && println("[Other tests skipped]\n")

# Use traditional Test.jl instead of TestItemRunner.jl for Aqua and JET
if (GROUP == "All") || (GROUP == "Code-Quality")
    println("Start running code quality tests...")

    Pkg.activate("core-test/code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using Test
    using QuantumToolbox
    using Aqua, JET

    include(joinpath(dirname(@__FILE__), "core-test", "code-quality", "code_quality.jl"))
end
