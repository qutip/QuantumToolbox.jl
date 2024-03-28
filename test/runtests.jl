using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if (GROUP == "All") || (GROUP == "Core")
    @time @safetestset "QuantumObjects" include("quantum_objects.jl")
    @time @safetestset "Time Evolution and partial trace" include("time_evolution_and_partial_trace.jl")
    @time @safetestset "Dynamical Fock Dimension mesolve" include("dynamical_fock_dimension_mesolve.jl")
    @time @safetestset "Dynamical Shifted Fock" include("dynamical-shifted-fock.jl")
    @time @safetestset "Generalized Master Equation" include("generalized_master_equation.jl")
    @time @safetestset "Eigenvalues and Operators" include("eigenvalues_and_operators.jl")
    @time @safetestset "Steady State" include("steady_state.jl")
    @time @safetestset "Entanglement" include("entanglement.jl")
    @time @safetestset "Negativity and Partial Transpose" include("negativity_and_partial_transpose.jl")
    @time @safetestset "Wigner" include("wigner.jl")
    @time @safetestset "Permutation" include("permutation.jl")
    @time @safetestset "Correlations and Spectrum" include("correlations_and_spectrum.jl")
    @time @safetestset "LowRankDynamics" include("low_rank_dynamics.jl")
end