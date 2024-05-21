using QuantumToolbox
using Documenter

DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, :(using QuantumToolbox); recursive = true)

const PAGES = Any[
    "Front Matter" => Any[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        # "Cite QuantumToolbox.jl" => "cite.md",
    ],
    "Users Guide" => Any[
        "Basic Operations on Quantum Objects" => Any[
            "Qobj.md",
            "Qobj_functions.md",
        ],
        "Manipulating States and Operators" => Any[
            "state_vectors.md",
            "density_matrices.md",
            "qubit_systems.md",
            "expectation_values.md",
            "superoperators.md",
        ],
        "Tensor Products and Partial Traces" => Any[
            "tensor.md",
            "partial_trace.md",
        ],
        "Time Evolution and Dynamics" => Any[
            "time_evolution/intro.md"
        ],
        "Solving for Steady-State Solutions" => Any[],
        "Two-time correlation functions" => Any[],
        "Extensions" => Any[
            "extensions/cuda.md",
        ],
    ],
    "Tutorials" => Any[
        "Time Evolution" => Any[
            "Low Rank Master Equation" => "tutorials/lowrank.md",
        ],
        "Miscellaneous Tutorials" => Any[
            "tutorials/logo.md",
        ],
    ],
    "Solver Benchmarks" => "benchmarks/benchmark_history.md",
    "API" => "api.md",
    # "Change Log" => "changelog.md",
]

makedocs(;
    modules = [QuantumToolbox],
    authors = "Alberto Mercurio",
    repo = Remotes.GitHub("albertomercurio", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    pages = PAGES,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://albertomercurio.github.io/QuantumToolbox.jl",
        edit_link = "main",
        assets = String[],
        mathengine = MathJax3(
            Dict(
                :loader => Dict("load" => ["[tex]/physics"]),
                :tex => Dict(
                    "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                    "tags" => "ams",
                    "packages" => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
    )
)

deploydocs(; repo = "github.com/albertomercurio/QuantumToolbox.jl", devbranch = "main")
