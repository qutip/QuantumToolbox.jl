#=
Some alias and overloading functions from QuantumToolboxCore
=#

function QuantumToolboxCore._add_library_deps!(lib::Val{:QuantumToolbox}, DEPpkgs::Vector{Module})
    QuantumToolboxCore._add_pkgs!(DEPpkgs, Module[SciMLBase, SciMLOperators, OrdinaryDiffEqCore, LinearSolve])
    return nothing
end

const settings = QuantumToolboxCore.settings
const versioninfo = QuantumToolboxCore.versioninfo
const about = QuantumToolboxCore.about
const cite = QuantumToolboxCore.cite
