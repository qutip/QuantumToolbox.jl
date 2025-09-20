#=
Command line output of information on QuantumToolbox, dependencies, and system information
=#

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
        "    Alberto Mercurio and Yi-Te Huang\n",
    )

    # print package information
    println(
        io,
        "Package information:\n",
        "====================================\n",
        "Julia              Ver. $(VERSION)\n",
        "QuantumToolbox     Ver. $(_get_pkg_version("QuantumToolbox"))\n",
        "SciMLOperators     Ver. $(_get_pkg_version("SciMLOperators"))\n",
        "LinearSolve        Ver. $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEqCore Ver. $(_get_pkg_version("OrdinaryDiffEqCore"))\n",
    )

    # print System information
    println(
        io,
        "System information:\n",
        "====================================\n",
        """OS       : $(OS_name) ($(Sys.MACHINE))\n""",
        """CPU      : $(length(cpu)) × $(cpu[1].model)\n""",
        """Memory   : $(round(Sys.total_memory() / 2 ^ 30, digits=3)) GB\n""",
        """WORD_SIZE: $(Sys.WORD_SIZE)\n""",
        """LIBM     : $(Base.libm_name)\n""",
        """LLVM     : libLLVM-$(Base.libllvm_version) ($(Sys.JIT), $(Sys.CPU_NAME))\n""",
        """BLAS     : $(basename(BLAS_info.libname)) ($(BLAS_info.interface))\n""",
        """Threads  : $(Threads.nthreads()) (on $(Sys.CPU_THREADS) virtual cores)\n""",
    )

    # print citation information
    println(
        io,
        "+---------------------------------------------------+\n",
        "| Please cite QuantumToolbox.jl in your publication |\n",
        "+---------------------------------------------------+\n",
        "For your convenience, a bibtex reference can be easily generated using `QuantumToolbox.cite()`.\n"
    )
    return nothing
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

@doc raw"""
    QuantumToolbox.cite(io::IO = stdout)

Command line output of citation information and bibtex generator for `QuantumToolbox.jl`.
"""
function cite(io::IO = stdout)
    citation = raw"""
    @article{QuantumToolbox-jl2025,
      doi = {},
      url = {},
      title = {{QuantumToolbox.jl}: {A}n efficient {J}ulia framework for simulating open quantum systems},
      author = {Mercurio, Alberto and Huang, Yi-Te and Cai, Li-Xun and Chen, Yueh-Nan and Savona, Vincenzo and Nori, Franco},
      journal = {{Quantum}},
      issn = {2521-327X},
      publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
      volume = {},
      pages = {},
      month = sep,
      year = {2025}
    }
    """

    return println(io, citation)
end
