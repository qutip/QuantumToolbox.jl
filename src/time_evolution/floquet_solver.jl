@doc raw"""
    struct FloquetProblem

Floquet counterpart to TimeEvolutionProblem for quantum systems with a periodic Hamitonian

# Fields:
- `prob::AbstractSciMLProblem`: The `ODEProblem` of the time evolution.
- `T::Real`: Period such that ``\hat{H}(T) = \hat{H}(0)``.
- `memoized_t::TT`: Times at which to pre-compute the system propagator and Floquet modes.
- `_tlist::TT`: Time list to use internally when computing the system propagator over a period T. Must contain all elements in `memoized_t`.
- `states_type::QuantumObjectType`: The type of the quantum states during the evolution (e.g., [`Ket`](@ref), [`Operator`](@ref), [`OperatorKet`](@ref), or [`SuperOperator`](@ref)).
- `dimensions::AbstractDimensions`: The dimensions of the Hilbert space.
- `kwargs::KWT`: Generic keyword arguments.

"""
struct FloquetProblem{
    ST<:QuantumObjectType,
    DT<:AbstractDimensions,
    PT<:AbstractSciMLProblem,
    TT<:AbstractVector,
    KWT,
} <: TimeEvolutionProblem
    prob::PT
    T<:Real
    memoized_t::TT
    _tlist::TT
    states_types::ST
    dimension::DT
    kwargs::KWT
end

"""
   struct FloquetBasis

Julia struct containing propagators, quasienergies, and Floquet states for a system with a ``T``-periodic Hamiltonain.

# Fields:
- `tlist::TT`: Times for which the system propagator and Floquet modes are precomputed. First and final elements are always ``0, T``.
- `H::AbstractQuantumObject`: T-periodic Hamiltonian.
- `T<:Real`: Hamiltonian period such that ``\hat{H}(T) = \hat{H}(0)``
- `U::Qobj`: System propagator at time ``T``
- `U_t::AbstractVector{Qobj}`: List of system propagators at times `tlist`
- `states::AbstractVector{<:AbstractQuantumObject}`: List of Floquet states at times `tlist`.
- `𝞊::AbstractVector{<:Real}`: Time-independent quasienergies
"""
struct FloquetBasis{
    TT<:AbstractVector{<:Real},
    TS<:AbstractVector,
    TU<:AbstractVector,
}
    tlist::TT
    H::AbstractQuantumObject
    T<:Real
    U::Qobj
    U_t::AbstractVector{Qobj}
    modes::AbstractVector{<:AbstractQuantumObject}
    𝞊::AbstractVector{<:Real}
end
