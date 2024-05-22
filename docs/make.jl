using QuantumToolbox
using Documenter

DocMeta.setdocmeta!(QuantumToolbox, :DocTestSetup, :(using QuantumToolbox); recursive = true)

makedocs(;
    modules = [QuantumToolbox],
    authors = ["Alberto Mercurio", "Luca Gravina", "Yi-Te Huang"],
    repo = Remotes.GitHub("qutip", "QuantumToolbox.jl"),
    sitename = "QuantumToolbox.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://qutip.github.io/QuantumToolbox.jl",
        edit_link = "main",
        assets = ["src/assets/favicon.ico"],
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
    ),
    pages = [
        "index.md",
        "api.md",
        "Tutorials" => [
            "Creation of the QuantulToolbox.jl logo" => "tutorials/logo.md",
            "Low Rank Master Equation" => "tutorials/lowrank.md",
        ],
        "Benchmarks" => ["Benchmark History" => "benchmarks/benchmark_history.md"],
    ],
)

deploydocs(; repo = "github.com/qutip/QuantumToolbox.jl", devbranch = "main")
