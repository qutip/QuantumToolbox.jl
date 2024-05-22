"""
    QuantumToolbox.versioninfo(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.about`](@ref).
"""
function versioninfo(io::IO = stdout)
    cpu = Sys.cpu_info()
    BLAS_info = BLAS.get_config().loaded_libs[1]
    Sys.iswindows() ? OS_name = "Windows" : Sys.isapple() ? OS_name = "macOS" : OS_name = Sys.KERNEL

    # print introduction
    println(
        io,
        "\n",
        " QuantumToolbox.jl: Quantum Toolbox in Julia\n",
        "≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡\n",
        "Copyright © QuTiP team 2022 and later.\n",
        "Current admin team:\n",
        "    Alberto Mercurio, Luca Gravina, Yi-Te Huang\n",
    )

    # print package informations
    println(
        io,
        "Package information:\n",
        "====================================\n",
        "QuantumToolbox  Ver. $(_get_pkg_version("QuantumToolbox"))\n",
        "LinearSolve     Ver. $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEq  Ver. $(_get_pkg_version("OrdinaryDiffEq"))\n",
    )

    # print System informations
    println(io, "System information:\n", "====================================\n", "Julia Version: $(VERSION)")
    println(io, """OS       : $(OS_name) ($(Sys.MACHINE))""")
    println(io, """CPU      : $(length(cpu)) × $(cpu[1].model)""")
    println(io, """Memory   : $(round(Sys.total_memory() / 2 ^ 30, digits=3)) GB""")
    println(io, """WORD_SIZE: $(Sys.WORD_SIZE)""")
    println(io, """LIBM     : $(Base.libm_name)""")
    println(io, """LLVM     : libLLVM-$(Base.libllvm_version) ($(Sys.JIT), $(Sys.CPU_NAME))""")
    println(io, """BLAS     : $(basename(BLAS_info.libname)) ($(BLAS_info.interface))""")
    println(io, """Threads  : $(Threads.nthreads()) (on $(Sys.CPU_THREADS) virtual cores)""")
    return print(io, "\n")
end

"""
    QuantumToolbox.about(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.versioninfo`](@ref).
"""
about(io::IO = stdout) = versioninfo(io)

function _get_pkg_version(pkg_name::String)
    D = Pkg.dependencies()
    for uuid in keys(D)
        if D[uuid].name == pkg_name
            return D[uuid].version
        end
    end
end
