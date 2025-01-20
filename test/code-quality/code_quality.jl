@testset "Code quality" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumToolbox; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        # JET.test_package(QuantumToolbox; target_defined_modules = true, ignore_missing_comparison = true)

        include("workload.jl")

        # Here we define some functins to exclude from the JET test

        # Related to the `check_error` function inside the `integrator` interface
        sci_ml_integrator_functions = (
            Base.modulesof!,
            Base.show,
            Base.show_at_namedtuple,
            Base.show_typealias,
            Base._show_type,
            Base.isvisible,
            Base.eltype,
        )

        # Related to FFTW.jl
        fftw_functions = (QuantumToolbox.unsafe_destroy_plan,)

        function_filter(@nospecialize f) = f ∉ (sci_ml_integrator_functions..., fftw_functions...)

        @test_opt function_filter = function_filter run_all_functions()
    end
end
