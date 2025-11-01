#! format: off
# turns off the julia formatting of this file

using QuantumToolbox
using Documenter
using DocumenterVitepress
using DocumenterCitations
using Changelog

# Load of packages required to compile the extension documentation
using CairoMakie

doctest_setup = quote
    using LinearAlgebra
    using SparseArrays
    using QuantumToolbox
end
DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, doctest_setup; recursive = true)

# some options for `makedocs`
const DRAFT = get(ENV, "DRAFT", false) == "true"  # `DRAFT   = true`  disables cell evaluation
const DOCTEST = get(ENV, "DOCTEST", true) == false # `DOCTEST = false` skips doc tests

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
        # "Key differences from QuTiP" => "getting_started/qutip_differences.md",
        "The Importance of Type-Stability" => "getting_started/type_stability.md",
        "Example: Create QuantumToolbox.jl Logo" => "getting_started/logo.md",
        "Cite QuantumToolbox.jl" => "getting_started/cite.md",
    ],
    "Users Guide" => [
        "Basic Operations on Quantum Objects" => [
            "Quantum Objects (Qobj)" => "users_guide/QuantumObject/QuantumObject.md",
            "Functions operating on Qobj" => "users_guide/QuantumObject/QuantumObject_functions.md",
        ],
        "Manipulating States and Operators" => "users_guide/states_and_operators.md",
        "Tensor Products and Partial Traces" => "users_guide/tensor.md",
        "Time Evolution and Dynamics" => [
            "Introduction" => "users_guide/time_evolution/intro.md",
            "Time Evolution Solutions" => "users_guide/time_evolution/solution.md",
            "Schrödinger Equation Solver" => "users_guide/time_evolution/sesolve.md",
            "Lindblad Master Equation Solver" => "users_guide/time_evolution/mesolve.md",
            "Monte Carlo Solver" => "users_guide/time_evolution/mcsolve.md",
            "Stochastic Solver" => "users_guide/time_evolution/stochastic.md",
            "Solving Problems with Time-dependent Hamiltonians" => "users_guide/time_evolution/time_dependent.md",
            "Bloch-Redfield master equation" => "users_guide/time_evolution/brmesolve.md",
        ],
        "Automatic Differentiation" => "users_guide/autodiff.md",
        "Intensive parallelization on a Cluster" => "users_guide/cluster.md",
        "Hierarchical Equations of Motion" => "users_guide/HEOM.md",
        "Solving for Steady-State Solutions" => "users_guide/steadystate.md",
        "Two-time correlation functions" => "users_guide/two_time_corr_func.md",
        "Plotting on the Bloch Sphere" => "users_guide/plotting_the_bloch_sphere.md",
        "QuantumToolbox Settings" => "users_guide/settings.md",
        "Extensions" => [
            "Extension for CUDA.jl" => "users_guide/extensions/cuda.md",
            "Extension for the Makie.jl ecosystem" => "users_guide/extensions/cairomakie.md",
        ],
    ],
    "Resources" => [
        "API" => "resources/api.md",
        "Bibliography" => "resources/bibliography.md",
        "ChangeLog" => "resources/changelog.md",
        "Contributing to QuantumToolbox.jl" => "resources/contributing.md",
        "Acknowledgements" => "resources/acknowledgements.md",
    ],
]

makedocs(;
    modules = [
        QuantumToolbox, 
        Base.get_extension(QuantumToolbox, :QuantumToolboxMakieExt),
    ],
    authors = "Alberto Mercurio and Yi-Te Huang",
    repo = Remotes.GitHub("qutip", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    pages = PAGES,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/qutip/QuantumToolbox.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    draft = DRAFT,
    doctest = DOCTEST,
    plugins = [bib],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/qutip/QuantumToolbox.jl",
    target = joinpath(@__DIR__, "build"),
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
