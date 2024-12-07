#! format: off
# turns off the julia formatting of this file

using QuantumToolbox
using Documenter
using DocumenterVitepress
using DocumenterCitations
using Changelog

DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, :(using QuantumToolbox); recursive = true)

const DRAFT = false  # set `true`  to disable cell evaluation
const DOCTEST = true # set `false` to skip doc tests

# generate bibliography
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "resources", "bibliography.bib"), 
    style=:authoryear,
)

# generate changelog
Changelog.generate(
    Changelog.Documenter(),
    joinpath(@__DIR__, "..", "CHANGELOG.md"),
    joinpath(@__DIR__, "src", "resources", "changelog.md");
    repo = "qutip/QuantumToolbox.jl",
)

const PAGES = [
    "Home" => "index.md",
    "Getting Started" => [
        "Brief Example" => "getting_started/brief_example.md",
        "Key differences from QuTiP" => "getting_started/qutip_differences.md",
        "The Importance of Type-Stability" => "getting_started/type_stability.md",
        # "Cite QuantumToolbox.jl" => "getting_started/cite.md",
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
        "Hierarchical Equations of Motion" => "users_guide/HEOM.md",
        "Solving for Steady-State Solutions" => "users_guide/steadystate.md",
        "Two-time correlation functions" => "users_guide/two_time_corr_func.md",
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
            "tutorials/cluster.md",
        ],
    ],
    "Resources" => [
        "API" => "resources/api.md",
        "Bibliography" => "resources/bibliography.md",
        "ChangeLog" => "resources/changelog.md",
        "Contributing to QuantumToolbox.jl" => "resources/contributing.md",
    ],
]

makedocs(;
    modules = [QuantumToolbox],
    authors = "Alberto Mercurio and Yi-Te Huang",
    repo = Remotes.GitHub("qutip", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    pages = PAGES,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/qutip/QuantumToolbox.jl",
    ),
    draft = DRAFT,
    doctest = DOCTEST,
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/qutip/QuantumToolbox.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
