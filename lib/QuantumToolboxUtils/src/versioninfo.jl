#=
Reusable version information helpers for QuantumToolbox libraries.
=#

# Registry of all loaded QuantumToolbox libraries, populated via __init__ in each library.
const _QT_LIBRARIES = Module[]

raw"""
    _register_qt_library!(m::Module)

Register a QuantumToolbox library module into the global registry. Each library should call this in its `__init__` function.
"""
function _register_qt_library!(m::Module)
    (m ∉ _QT_LIBRARIES) && pushfirst!(_QT_LIBRARIES, m) # use pushfirst! so that main API libraries are at the front of the registry (for better display order in versioninfo)
    return nothing
end

raw"""
    _add_library_deps!(lib::Val, DEPpkgs::Vector{Module})

Add new dependencies to `DEPpkgs` for the specified library. Extend this function in each QuantumToolbox library to declare its external (non-QuantumToolbox) dependencies.
"""
_add_library_deps!(lib::Val{T}, DEPpkgs::Vector{Module}) where {T} = throw(ArgumentError("Unknown QuantumToolbox library : $T"))
function _add_library_deps!(lib::Val{:QuantumToolboxUtils}, DEPpkgs::Vector{Module})
    _add_pkgs!(DEPpkgs, Module[SciMLOperators])
    return nothing
end

# add pkgs2 into pkgs1 but ensure the uniqueness
function _add_pkgs!(pkgs1::Vector{Module}, pkgs2::Vector{Module})
    for pkg in unique(pkgs2)
        (pkg ∉ pkgs1) && push!(pkgs1, pkg)
    end
    return nothing
end

raw"""
    QuantumToolboxUtils.pkginfo(io::IO=stdout; pkgs::Vector{Module} = Module[])

Command line output of version numbers for given vector of packages: `pkgs`.
"""
function pkginfo(io::IO = stdout; pkgs::Vector{Module} = Module[], split_after::Union{Nothing, Int} = nothing)
    pkg_ver_list = map(pkgversion, pkgs)
    maxLen = max(5, maximum(length ∘ string, pkgs)) # maximum string length of package names (5 refer to "Julia")

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
    QuantumToolboxUtils.sysinfo(io::IO=stdout)

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

    DEPpkgs = _gen_dep_pkg_list()
    pkginfo(io; pkgs = vcat(_QT_LIBRARIES, DEPpkgs), split_after = length(_QT_LIBRARIES))
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

function _gen_dep_pkg_list()
    DEPpkgs = Module[]
    for lib in _QT_LIBRARIES
        _add_library_deps!(Val(nameof(lib)), DEPpkgs)
    end
    return DEPpkgs
end

@doc raw"""
    QuantumToolbox.versioninfo(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.about`](@ref QuantumToolboxUtils.about).
"""
versioninfo(io::IO = stdout) = _print_versioninfo(io)

@doc raw"""
    QuantumToolbox.about(io::IO=stdout)

Command line output of information on QuantumToolbox, dependencies, and system information, same as [`QuantumToolbox.versioninfo`](@ref QuantumToolboxUtils.versioninfo).
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
