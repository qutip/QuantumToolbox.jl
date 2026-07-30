#=
Reusable version information helpers for QuantumToolbox libraries.
=#

# Registry of all loaded QuantumToolbox libraries, populated via __init__ in each library.
const QT_LIBRARIES = Module[]
const DEP_PKGS = Module[]

raw"""
    QuantumToolboxCore.pkginfo(io::IO=stdout; pkgs::Vector{Module} = Module[])

Command line output of version numbers for given vector of packages: `pkgs`.
"""
function pkginfo(io::IO = stdout; pkgs::Vector{Module} = Module[], split_after::Union{Nothing, Int} = nothing)
    pkg_ver_list = map(pkgversion, pkgs)

    # maximum string length of package names (5 refer to "Julia")
    maxLen = isempty(pkgs) ? max(5, maximum(length ∘ string, pkgs)) : 5

    separation_line = "------------------------------------"
    print(
        io,
        "Package information:\n",
        "====================================\n",
    )
    println(io, rpad("Julia", maxLen, " "), " Ver. ", VERSION) # print Julia version first
    println(io, separation_line)
    for (idx, (pkg, pkg_ver)) in enumerate(zip(pkgs, pkg_ver_list))
        println(io, rpad(pkg, maxLen, " "), " Ver. ", pkg_ver)
        !isnothing(split_after) && (idx == split_after) && (idx < length(pkgs)) && println(io, separation_line)
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
        "====================================\n",
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

    pkginfo(io; pkgs = vcat(QT_LIBRARIES, DEP_PKGS), split_after = length(QT_LIBRARIES))
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
