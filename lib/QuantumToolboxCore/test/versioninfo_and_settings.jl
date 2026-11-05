using Test
using QuantumToolboxCore

@testset "versioninfo and about" begin
    # versioninfo
    io_version = IOBuffer()
    QuantumToolboxCore.versioninfo(io_version)
    version_output = String(take!(io_version))
    @test occursin("QuantumToolbox.jl: Quantum Toolbox in Julia", version_output)
    @test occursin("Package information:", version_output)
    @test occursin("QuantumToolboxCore", version_output)
    @test occursin("System information:", version_output)
    @test occursin("Please cite QuantumToolbox.jl in your publication", version_output)

    # about
    io_about = IOBuffer()
    QuantumToolboxCore.about(io_about)
    about_output = String(take!(io_about))
    @test version_output == about_output
end

@testset "Settings" begin
    settings = QuantumToolboxCore.settings

    io = IOBuffer()
    show(io, settings)
    out = String(take!(io))
    @test occursin("QuantumToolbox.jl Settings", out)
    for name in propertynames(settings)
        @test occursin(string(name), out)
    end

    # these settings change a type, so reading them must stay type-stable
    get_eltype() = settings.default_eltype
    get_kwargs() = settings.progressmeter_kwargs
    @test @inferred(get_eltype()) === ComplexF64
    @test @inferred(get_kwargs()) == (showspeed = true, printed = true)

    try
        @testset "default_eltype" begin
            @test eltype(@inferred(destroy(3))) === ComplexF64

            settings.default_eltype = ComplexF32
            @test eltype(@inferred(destroy(3))) === ComplexF32
            @test eltype(@inferred(fock(3, 1))) === ComplexF32
            @test eltype(@inferred(sigmax())) === ComplexF32
            @test eltype(@inferred(rand_ket(3))) === ComplexF32
            @test eltype(@test_logs((:warn, r"higher precision"), destroy(ComplexF64, 3))) === ComplexF64 # an explicit type still wins

            settings.default_eltype = Float32
            @test eltype(destroy(3)) === Float32
            @test eltype(rand_ket(3)) === ComplexF32 # only defined for complex types
            @test eltype(rand_dm(3)) === ComplexF32
            @test eltype(qft(2)) === ComplexF32
            @test eltype(QuantumToolboxCore.momentum(3)) === ComplexF32

            for T in (Int, Complex{Int}, AbstractFloat, Complex, String)
                @test_throws ArgumentError settings.default_eltype = T
            end
            @test settings.default_eltype === Float32

            # changing it inside a function: the function cannot see the new value, so it throws instead
            change_and_use() = (settings.default_eltype = ComplexF16; destroy(3))
            @test_throws ErrorException change_and_use()
            @test eltype(destroy(3)) === ComplexF16
        end

        @testset "default_eltype_warning" begin
            settings.default_eltype = ComplexF32
            @test_logs (:warn, r"higher precision") QuantumObject(rand(ComplexF64, 2))
            @test_logs (:warn, r"higher precision") 1.0 * destroy(3)
            @test_logs QuantumObject(rand(ComplexF32, 2))
            @test_logs QuantumObject(rand(Float32, 2))
            @test_logs QuantumObject(rand(Int, 2))
            settings.default_eltype_warning = false
            @test_logs QuantumObject(rand(ComplexF64, 2))
            A = QuantumObject(rand(ComplexF64, 2, 2))
            settings.default_eltype_warning = true
            @test_logs (:warn, r"higher precision") QobjEvo((A, (p, t) -> sin(t)))

            # with defaults that `Float64` numbers cannot promote, a higher precision is always chosen explicitly
            settings.default_eltype = ComplexF64
            @test_logs QuantumObject(rand(Complex{BigFloat}, 2))
        end

        @testset "progressmeter_kwargs" begin
            settings.progressmeter_kwargs = (showspeed = false, output = devnull)
            @test @inferred(get_kwargs()) == (showspeed = false, output = devnull)
            @test_throws ArgumentError settings.progressmeter_kwargs = 1

            # deprecated name
            @test (@test_logs (:warn, r"deprecated") settings.ProgressMeterKWARGS) == (showspeed = false, output = devnull)
            @test_logs (:warn, r"deprecated") settings.ProgressMeterKWARGS = (printed = true,)
            @test settings.progressmeter_kwargs == (printed = true,)
        end

        @testset "tidyup_tol and auto_tidyup" begin
            # these are plain fields: changing them inside a function works
            change_tol() = (settings.tidyup_tol = 1.0e-3; settings.tidyup_tol)
            @test change_tol() == 1.0e-3
            @test_throws MethodError settings.tidyup_tol = "1e-3"
        end
    finally
        settings.tidyup_tol = 1.0e-14
        settings.auto_tidyup = true
        settings.default_eltype = ComplexF64
        settings.default_eltype_warning = true
        settings.progressmeter_kwargs = (showspeed = true, printed = true)
    end
    @test eltype(destroy(3)) === ComplexF64
end
