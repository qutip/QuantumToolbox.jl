using Test
using QuantumToolbox
using Aqua
using JET

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

if (GROUP == "All") || (GROUP == "Code Quality")
    using Aqua
    using JET
    include(joinpath(testdir, "aqua.jl"))
    include(joinpath(testdir, "jet.jl"))
end