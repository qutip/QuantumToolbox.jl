#=
Some alias and overloading functions from QuantumToolboxUtils
=#

function QuantumToolboxUtils._add_library_deps!(lib::Val{:QuantumToolbox}, DEPpkgs::Vector{Module})
    QuantumToolboxUtils._add_pkgs!(DEPpkgs, Module[SciMLBase, SciMLOperators, OrdinaryDiffEqCore, LinearSolve])
    return nothing
end

const settings = QuantumToolboxUtils.settings
const versioninfo = QuantumToolboxUtils.versioninfo
const about = QuantumToolboxUtils.about
const cite = QuantumToolboxUtils.cite
