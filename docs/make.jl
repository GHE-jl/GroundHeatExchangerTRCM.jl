using ThermalResistanceCapacitanceModel
using Documenter

# Make `using ThermalResistanceCapacitanceModel` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    ThermalResistanceCapacitanceModel,
    :DocTestSetup,
    :(using ThermalResistanceCapacitanceModel);
    recursive = true,
)

makedocs(;
    modules = [ThermalResistanceCapacitanceModel],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "ThermalResistanceCapacitanceModel.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/ThermalResistanceCapacitanceModel.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => [
            "Cross-section networks" => "theory/cross_section.md",
            "Response model" => "theory/response_model.md",
            "Short and long term" => "theory/hybrid.md",
            "Borehole-to-borehole responses" => "theory/pair_response.md",
        ],
        "Design notes" => "notes.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/ThermalResistanceCapacitanceModel.jl",
    devbranch = "main",
)
