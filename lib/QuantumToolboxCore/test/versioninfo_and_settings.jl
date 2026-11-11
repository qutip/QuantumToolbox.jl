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
    @test occursin("default_eltype", sprint(show, settings))

    # `default_eltype` changes a type, so reading it must be type-stable
    get_eltype() = settings.default_eltype
    @test @inferred(get_eltype()) === ComplexF64

    settings.default_eltype = ComplexF32
    @test eltype(@inferred(destroy(3))) === ComplexF32
    settings.default_eltype = Float32
    @test eltype(destroy(3)) === Float32
    @test eltype(@test_logs((:warn, r"more general"), rand_ket(3))) === ComplexF32 # only defined for complex types
    @test_throws ArgumentError settings.default_eltype = Int

    # the function cannot see the value it has just changed, so it throws instead of using the previous one
    change_and_use() = (settings.default_eltype = ComplexF32; destroy(3))
    @test_throws ErrorException change_and_use()

    # warn when the element type is more general than the default one
    @test_logs (:warn, r"more general") 0.5 * destroy(3)
    @test_logs (:warn, r"more general") QuantumObject(rand(Complex{BigFloat}, 2))
    @test_logs QuantumObject(rand(Float16, 2))
    settings.show_eltype_warning = false
    @test_logs 0.5 * destroy(3)

    settings.progressmeter_kwargs = (showspeed = false, output = devnull)
    @test settings.progressmeter_kwargs == (showspeed = false, output = devnull)
    @test (@test_logs (:warn, r"deprecated") settings.ProgressMeterKWARGS) == settings.progressmeter_kwargs

    # restore the defaults, so that the other tests are not affected
    settings.default_eltype = ComplexF64
    settings.show_eltype_warning = true
    settings.progressmeter_kwargs = (showspeed = true, printed = true)
end
