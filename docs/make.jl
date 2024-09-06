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
        "Examples" => [
            joinpath("examples", "setup.md"),
            joinpath("examples", "gwas.md"),
            joinpath("examples", "phewas.md"),
            joinpath("examples", "interactions.md"),
        ],
        "User Guide" => [
            "overview.md",
            "all_workflows_parameters.md",
            "runtime_considerations.md",
            "The Discovery Workflow" => [
                joinpath("targene", "overview.md"),
                joinpath("targene", "data_sources.md"),
                joinpath("targene", "confounding_adjustment.md"),
                joinpath("targene", "study_designs.md"),
                joinpath("targene", "tmle.md"),
                joinpath("targene", "sieve_variance.md"),
                joinpath("targene", "miscellaneous.md"),
                joinpath("targene", "outputs.md"),
            ],
            "The Simulation Workflows (Experimental)" => [
                joinpath("simulations", "overview.md"),
                joinpath("simulations", "null_simulation.md"),
                joinpath("simulations", "realistic_simulation.md"),
            ],
            ],
        "Secondary Workflows" => [
                joinpath("secondary_workflows", "pca.md"),
                joinpath("secondary_workflows", "make_dataset.md")
            ],
        "Developper Guide" => [
            joinpath("developer_guide", "project_organization.md"),
            joinpath("developer_guide", "contribution_guide.md")
            ],
    ],
)

deploydocs(;
    repo="github.com/TARGENE/targene-pipeline",
    devbranch="main",
    push_preview=true
)
