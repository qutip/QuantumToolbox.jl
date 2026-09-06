#=
Reusable version information helpers for QuantumToolbox libraries.
=#

# Registry of all loaded QuantumToolbox libraries, populated via __init__ in each library.
const QT_LIBRARIES = Module[]
const DEP_PKGS = Module[]
const EXT_PKGS = Module[]

# separation lines
const SEPARATION_LINE_LENGTH = 36
const SINGLE_SEPARATION_LINE = repeat("-", SEPARATION_LINE_LENGTH) * "\n"
const DOUBLE_SEPARATION_LINE = repeat("=", SEPARATION_LINE_LENGTH) * "\n"

raw"""
    QuantumToolboxCore.pkginfo(io::IO=stdout)

Command line output of version numbers for:

- QuantumToolbox libraries
- Dependencies
- Triggered extensions
"""
function pkginfo(io::IO = stdout)

    all_pkgs = vcat(QT_LIBRARIES, DEP_PKGS, EXT_PKGS)
    all_pkgs_ver = map(pkgversion, all_pkgs)

    # maximum string length of package names (5 refer to "Julia")
    maxLen = max(5, maximum(length ∘ string, all_pkgs))

    print(
        io,
        "Package information:\n",
        DOUBLE_SEPARATION_LINE,
    )
    println(io, rpad("Julia", maxLen, " "), " Ver. ", VERSION) # print Julia version first

    idx = 1  # index for all_pkgs_ver iteration

    # QuantumToolbox libraries
    for pkg in QT_LIBRARIES
        println(io, rpad(pkg, maxLen, " "), " Ver. ", all_pkgs_ver[idx])
        idx += 1
    end

    # dependencies
    if !isempty(DEP_PKGS)
        println(io, SINGLE_SEPARATION_LINE, "dependencies:")
        for pkg in DEP_PKGS
            println(io, rpad(pkg, maxLen, " "), " Ver. ", all_pkgs_ver[idx])
            idx += 1
        end
    end

    # triggered extensions
    if !isempty(EXT_PKGS)
        println(io, SINGLE_SEPARATION_LINE, "triggered extensions for:")
        for pkg in EXT_PKGS
            println(io, rpad(pkg, maxLen, " "), " Ver. ", all_pkgs_ver[idx])
            idx += 1
        end
    end

    print(io, "\n")
    return nothing
end

raw"""
    QuantumToolboxCore.sysinfo(io::IO=stdout)

Command line output of system information.
"""
function sysinfo(io::IO = stdout)
    cpu = Sys.cpu_info()
    BLAS_info = LinearAlgebra.BLAS.get_config().loaded_libs[1]
    Sys.iswindows() ? OS_name = "Windows" : Sys.isapple() ? OS_name = "macOS" : OS_name = Sys.KERNEL

    println(
        io,
        "System information:\n",
        DOUBLE_SEPARATION_LINE,
        """OS       : $(OS_name) ($(Sys.MACHINE))\n""",
        """CPU      : $(length(cpu)) × $(cpu[1].model)\n""",
        """Memory   : $(round(Sys.total_memory() / 2^30, digits = 3)) GB\n""",
        """WORD_SIZE: $(Sys.WORD_SIZE)\n""",
        """LIBM     : $(Base.libm_name)\n""",
        """LLVM     : libLLVM-$(Base.libllvm_version) ($(Sys.JIT), $(Sys.CPU_NAME))\n""",
        """BLAS     : $(basename(BLAS_info.libname)) ($(BLAS_info.interface))\n""",
        """Threads  : $(Threads.nthreads()) (on $(Sys.CPU_THREADS) virtual cores)\n""",
    )
    return nothing
end

function _print_versioninfo(io::IO = stdout)
    println(
        io,
        "\n",
        " QuantumToolbox.jl: Quantum Toolbox in Julia\n",
        "≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡\n",
        "Copyright © QuTiP team 2022 and later.\n",
        "Current admin team:\n",
        "    Alberto Mercurio and Yi-Te Huang\n",
    )

    pkginfo(io)
    sysinfo(io)

    println(
        io,
        "+---------------------------------------------------+\n",
        "| Please cite QuantumToolbox.jl in your publication |\n",
        "+---------------------------------------------------+\n",
        "For your convenience, a bibtex reference can be easily generated using `QuantumToolbox.cite()`.\n",
    )
    return nothing
end

@doc raw"""
    QuantumToolbox.versioninfo(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.about`](@ref QuantumToolboxCore.about).
"""
versioninfo(io::IO = stdout) = _print_versioninfo(io)

@doc raw"""
    QuantumToolbox.about(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.versioninfo`](@ref QuantumToolboxCore.versioninfo).
"""
about(io::IO = stdout) = versioninfo(io)

@doc raw"""
    QuantumToolbox.cite(io::IO = stdout)

Command line output of citation information and bibtex generator for `QuantumToolbox.jl`.
"""
function cite(io::IO = stdout)
    citation = raw"""
    @article{QuantumToolbox.jl2025,
      title = {Quantum{T}oolbox.jl: {A}n efficient {J}ulia framework for simulating open quantum systems},
      author = {Mercurio, Alberto and Huang, Yi-Te and Cai, Li-Xun and Chen, Yueh-Nan and Savona, Vincenzo and Nori, Franco},
      journal = {{Quantum}},
      issn = {2521-327X},
      publisher = {{Verein zur F{\\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
      volume = {9},
      pages = {1866},
      month = sep,
      year = {2025},
      doi = {10.22331/q-2025-09-29-1866},
      url = {https://doi.org/10.22331/q-2025-09-29-1866}
    }
    """
    return print(io, citation)
end
