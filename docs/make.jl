#! format: off
# turns off the julia formatting of this file

using QuantumToolbox
using Documenter
using DocumenterVitepress
using DocumenterCitations

DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, :(using QuantumToolbox); recursive = true)

const DRAFT = false # set `true` to disable cell evaluation

bib = CitationBibliography(joinpath(@__DIR__, "src", "bibliography.bib"), style=:authoryear)

# const MathEngine = MathJax3(
#     Dict(
#         :loader => Dict("load" => ["[tex]/physics"]),
#         :tex => Dict(
#             "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
#             "tags" => "ams",
#             "packages" => ["base", "ams", "autoload", "physics"],
#         ),
#     )
# )

const PAGES = [
    "Getting Started" => [
        "Introduction" => "index.md",
        "Key differences from QuTiP" => "qutip_differences.md",
        "The Importance of Type-Stability" => "type_stability.md",
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
            "Introduction" => "users_guide/time_evolution/intro.md",
            "Time Evolution Solutions" => "users_guide/time_evolution/solution.md",
            "Schrödinger Equation Solver" => "users_guide/time_evolution/sesolve.md",
            "Lindblad Master Equation Solver" => "users_guide/time_evolution/mesolve.md",
            "Monte-Carlo Solver" => "users_guide/time_evolution/mcsolve.md",
            "Stochastic Solver" => "users_guide/time_evolution/stochastic.md",
            "Solving Problems with Time-dependent Hamiltonians" => "users_guide/time_evolution/time_dependent.md",
        ],
        "Solving for Steady-State Solutions" => "users_guide/steadystate.md",
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
    "Bibliography" => "bibliography.md",
    # "Change Log" => "changelog.md",
]

makedocs(;
    modules = [QuantumToolbox],
    authors = "Alberto Mercurio and Yi-Te Huang",
    repo = Remotes.GitHub("qutip", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    pages = PAGES,
    # format = Documenter.HTML(;
    #     prettyurls = get(ENV, "CI", "false") == "true",
    #     canonical = "https://qutip.github.io/QuantumToolbox.jl",
    #     edit_link = "main",
    #     assets = ["assets/favicon.ico"],
    #     mathengine = MathEngine,
    #     size_threshold_ignore = ["api.md"],
    # ),
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://qutip.github.io/QuantumToolbox.jl",
        # deploy_url = "https://qutip.org/QuantumToolbox.jl/",
    ),
    draft = DRAFT,
    plugins = [bib],
)

# deploydocs(; repo = "github.com/qutip/QuantumToolbox.jl", devbranch = "main")

deploydocs(;
    repo = "github.com/qutip/QuantumToolbox.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
