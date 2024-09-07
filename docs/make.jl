#! format: off
# turns off the julia formatting of this file

using QuantumToolbox
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, :(using QuantumToolbox); recursive = true)

const MathEngine = MathJax3(
    Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
            "tags" => "ams",
            "packages" => ["base", "ams", "autoload", "physics"],
        ),
    )
)

const PAGES = [
    "Getting Started" => [
        "Introduction" => "index.md",
        "Key differences from QuTiP" => "qutip_differences.md",
        # "Cite QuantumToolbox.jl" => "cite.md",
    ],
    "Users Guide" => [
        "Basic Operations on Quantum Objects" => [
            "users_guide/QuantumObject/QuantumObject.md",
            "users_guide/QuantumObject/QuantumObject_functions.md",
        ],
        "Manipulating States and Operators" => "users_guide/states_and_operators.md",
        "Tensor Products and Partial Traces" => "users_guide/tensor.md",
        "Time Evolution and Dynamics" => [
            # "Introduction" => "users_guide/time_evolution/intro.md",
        ],
        "Solving for Steady-State Solutions" => [],
        "Symmetries" => [],
        "Two-time correlation functions" => [],
        "Extensions" => [
            "users_guide/extensions/cuda.md",
        ],
    ],
    "Tutorials" => [
        "Time Evolution" => [
            "Low Rank Master Equation" => "tutorials/lowrank.md",
        ],
        "Miscellaneous Tutorials" => [
            "tutorials/logo.md",
        ],
    ],
    "API" => "api.md",
    # "Change Log" => "changelog.md",
]

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules = [QuantumToolbox],
    authors = "Alberto Mercurio, Luca Gravina and Yi-Te Huang",
    repo = Remotes.GitHub("qutip", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    pages = PAGES,
    plugins=[bib],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://qutip.github.io/QuantumToolbox.jl",
        edit_link = "main",
        assets = ["assets/favicon.ico"],
        mathengine = MathEngine,
        size_threshold_ignore = ["api.md"],
    )
)

deploydocs(; repo = "github.com/qutip/QuantumToolbox.jl", devbranch = "main")
