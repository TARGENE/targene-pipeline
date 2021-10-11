using UKBBEpistasisPipeline
using Documenter

DocMeta.setdocmeta!(UKBBEpistasisPipeline, :DocTestSetup, :(using UKBBEpistasisPipeline); recursive=true)

makedocs(;
    modules=[UKBBEpistasisPipeline],
    authors="Olivier Labayle",
    repo="https://github.com/olivierlabayle/UKBBEpistasisPipeline.jl/blob/{commit}{path}#{line}",
    sitename="UKBBEpistasisPipeline.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olivierlabayle.github.io/UKBBEpistasisPipeline.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olivierlabayle/UKBBEpistasisPipeline.jl",
)
