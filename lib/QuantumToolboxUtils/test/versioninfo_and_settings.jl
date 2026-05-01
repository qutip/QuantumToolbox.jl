@testitem "versioninfo and about" begin
    # test error handling
    io_wrong = IOBuffer()
    push!(QuantumToolboxUtils._QT_LIBRARIES, Base)
    @test_throws ArgumentError QuantumToolboxUtils.versioninfo(io_wrong)
    deleteat!(QuantumToolboxUtils._QT_LIBRARIES, findall(x -> x == Base, QuantumToolboxUtils._QT_LIBRARIES))

    # versioninfo
    io_version = IOBuffer()
    QuantumToolboxUtils.versioninfo(io_version)
    version_output = String(take!(io_version))
    @test occursin("QuantumToolbox.jl: Quantum Toolbox in Julia", version_output)
    @test occursin("Package information:", version_output)
    @test occursin("QuantumToolboxUtils", version_output)
    @test occursin("System information:", version_output)
    @test occursin("Please cite QuantumToolbox.jl in your publication", version_output)

    # about
    io_about = IOBuffer()
    QuantumToolboxUtils.about(io_about)
    about_output = String(take!(io_about))
    @test version_output == about_output
end

@testitem "Settings" begin
    io = IOBuffer()
    show(io, QuantumToolboxUtils.settings)
    out = String(take!(io))
    @test occursin("QuantumToolbox.jl Settings", out)
    @test occursin("tidyup_tol", out)
    @test occursin("auto_tidyup", out)
    @test occursin("ProgressMeterKWARGS", out)
end
