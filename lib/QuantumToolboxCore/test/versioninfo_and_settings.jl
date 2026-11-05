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

    # `default_eltype` and `progressmeter_kwargs` change a type, so reading them must be type-stable
    get_eltype() = settings.default_eltype
    get_kwargs() = settings.progressmeter_kwargs
    @test @inferred(get_eltype()) === ComplexF64
    @test @inferred(get_kwargs()) == (showspeed = true, printed = true)

    # they are saved as preferences of the active project, and used after restarting Julia
    active_project = Base.ACTIVE_PROJECT[]
    mktempdir() do dir
        Base.ACTIVE_PROJECT[] = joinpath(dir, "Project.toml")
        try
            @test_logs (:info, r"Restart Julia") settings.default_eltype = ComplexF32
            @test_logs (:info, r"Restart Julia") settings.progressmeter_kwargs = (showspeed = false, output = devnull)
            prefs = read(joinpath(dir, "LocalPreferences.toml"), String)
            @test occursin("default_eltype = \"ComplexF32\"", prefs)
            @test occursin("progressmeter_kwargs = \"(showspeed = false, output = Base.DevNull())\"", prefs)
            @test settings.default_eltype === ComplexF64 # still the one in use
        finally
            Base.ACTIVE_PROJECT[] = active_project
        end
    end
    @test_throws ArgumentError settings.default_eltype = Int
    @test_throws ArgumentError settings.progressmeter_kwargs = (output = stdout,) # it cannot be read back
    @test (@test_logs (:warn, r"deprecated") settings.ProgressMeterKWARGS) == settings.progressmeter_kwargs

    # warn when the element type is more general than the default one
    @test_logs (:warn, r"more general") QuantumObject(rand(Complex{BigFloat}, 2))
    @test_logs QuantumObject(rand(ComplexF32, 2))
    settings.default_eltype_warning = false
    @test_logs QuantumObject(rand(Complex{BigFloat}, 2))
    settings.default_eltype_warning = true
end
