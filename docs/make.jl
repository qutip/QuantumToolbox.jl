using QuPhys
using Documenter

DocMeta.setdocmeta!(QuPhys, :DocTestSetup, :(using QuPhys); recursive=true)

makedocs(;
    modules=[QuPhys],
    authors="Alberto Mercurio",
    repo="https://github.com/albertomercurio/QuPhys.jl/blob/{commit}{path}#{line}",
    sitename="QuPhys.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://albertomercurio.github.io/QuPhys.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/albertomercurio/QuPhys.jl",
    devbranch="main",
)
