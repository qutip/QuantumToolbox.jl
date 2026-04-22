#=
Overloading functions from QuantumToolboxUtils
=#

function QuantumToolboxUtils._add_library_deps!(lib::Val{:QuantumToolbox}, DEPpkgs::Vector{Module})
    QuantumToolboxUtils._add_pkgs!(DEPpkgs, Module[SciMLBase, SciMLOperators, OrdinaryDiffEqCore, LinearSolve])
    return nothing
end

@doc raw"""
    QuantumToolbox.versioninfo(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.about`](@ref).
"""
versioninfo(io::IO = stdout) = QuantumToolboxUtils._print_versioninfo(io)

@doc raw"""
    QuantumToolbox.about(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.versioninfo`](@ref).
"""
about(io::IO = stdout) = versioninfo(io)

@doc raw"""
    QuantumToolbox.cite(io::IO = stdout)

Command line output of citation information and bibtex generator for `QuantumToolbox.jl`.
"""
cite(io::IO = stdout) = QuantumToolboxUtils.cite(io)
