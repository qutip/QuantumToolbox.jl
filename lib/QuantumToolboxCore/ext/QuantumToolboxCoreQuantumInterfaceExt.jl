module QuantumToolboxCoreQuantumInterfaceExt

#=
Some functions exported by QuantumToolbox have the same names as the ones owned by QuantumInterface.jl, which are
exported by the packages built on it (e.g., QuantumOptics.jl, SecondQuantizedAlgebra.jl, QuantumClifford.jl, and Gabs.jl).
This causes a name conflict when these packages are loaded together with QuantumToolbox, which the user can resolve by
explicitly importing the functions from one of them (e.g., `using SecondQuantizedAlgebra: tensor, ⊗`).

This extension does NOT change the functions exported by QuantumToolbox. It only defines the QuantumInterface.jl
functions for the QuantumToolbox types, so that they still work on QuantumObjects after such an explicit import.

Note that `tensor` and `⊗` of QuantumInterface.jl (v0.4.4 and above) are the ones owned by TensorCore.jl.

The following conflicting functions are not extended on purpose:
- `basis`, `displace`, `squeeze`: the QuantumToolbox methods only take numbers as arguments, so it would be type piracy.
- `ptrace`: QuantumToolbox takes the subsystems to keep, while QuantumInterface.jl takes the subsystems to trace out.
=#

using QuantumToolboxCore
import LinearAlgebra: kron
import QuantumInterface

QuantumInterface.tensor(A::AbstractQuantumObject, B::AbstractQuantumObject...) = kron(A, B...)
QuantumInterface.tensor(A::Vector{<:AbstractQuantumObject}) = kron(A)

QuantumInterface.expect(O::AbstractQuantumObject, ψ) = expect(O, ψ)
QuantumInterface.expect(O::AbstractVector{<:AbstractQuantumObject}, ψ) = expect(O, ψ)

QuantumInterface.variance(O::QuantumObject, ψ) = variance(O, ψ)

QuantumInterface.fidelity(ρ::QuantumObject, σ::QuantumObject) = fidelity(ρ, σ)

QuantumInterface.entropy_vn(ρ::QuantumObject; kwargs...) = entropy_vn(ρ; kwargs...)

# `QuantumInterface.commutator` is only available in QuantumInterface.jl v0.4.2 and above
if isdefined(QuantumInterface, :commutator)
    QuantumInterface.commutator(A::QuantumObject, B::QuantumObject; kwargs...) = commutator(A, B; kwargs...)
end

end
