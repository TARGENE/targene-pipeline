using Documenter

# DocMeta.setdocmeta!(TarGene, :DocTestSetup, :(using TarGene); recursive=true)

makedocs(;
    # modules=[TarGene],
    authors="Olivier Labayle",
    repo="https://github.com/TARGENE/targene-pipeline/blob/{commit}{path}#{line}",
    sitename="TarGene",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TARGENE.github.io/targene-pipeline",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => [
            "overview.md",
            "data_sources.md",
            "confounding_adjustment.md",
            "study_designs.md",
            "tmle.md",
            "sieve_variance.md",
            "miscellaneous.md",
            "negative_control.md"
        ],
        "Miscellaneous" => [
            "nextflow_params.md",
            "runtime_considerations.md"
        ],
        "Developper Guide" => [
            "project_organization.md",
            "contribution_guide.md"
            ],
        "Associated Softwares" => "associated_softwares.md",
        "Related Publications" => "publications.md"
    ],
)

deploydocs(;
    repo="github.com/TARGENE/targene-pipeline",
    devbranch="main",
    push_preview=true
)
