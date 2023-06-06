using QuPhys
using Documenter

ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

DocMeta.setdocmeta!(QuPhys, :DocTestSetup, :(using QuPhys); recursive=true)

makedocs(;
    modules=[QuPhys],
    authors="Alberto Mercurio",
    repo="https://github.com/albertomercurio/QuPhys.jl/blob/{commit}{path}#{line}",
    sitename="QuPhys.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://albertomercurio.github.io/QuPhys.jl",
        edit_link="main",
        assets=String[],
        mathengine = MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/physics"]),
            :tex => Dict(
                "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
                "tags" => "ams",
                "packages" => ["base", "ams", "autoload", "physics"],
            ),
            )),
    ),
    pages=[
        "index.md",
        "api.md",
        "Benchmark Results" => "../../benchmark/results.md",
    ],
)

deploydocs(;
    repo="github.com/albertomercurio/QuPhys.jl",
    devbranch="main",
)
