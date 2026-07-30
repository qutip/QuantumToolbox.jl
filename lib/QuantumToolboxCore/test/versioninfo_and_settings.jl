@testitem "versioninfo and about" begin
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

@testitem "Settings" begin
    io = IOBuffer()
    show(io, QuantumToolboxCore.settings)
    out = String(take!(io))
    @test occursin("QuantumToolbox.jl Settings", out)
    @test occursin("tidyup_tol", out)
    @test occursin("auto_tidyup", out)
    @test occursin("ProgressMeterKWARGS", out)
end
