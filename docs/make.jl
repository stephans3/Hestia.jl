using Hestia
using Documenter

DocMeta.setdocmeta!(Hestia, :DocTestSetup, :(using Hestia); recursive=true)

makedocs(;
    modules=[Hestia],
    authors="Stephan Scholz <stephan.scholz@rwu.de> and contributors",
    repo="https://github.com/stephans3/Hestia.jl/blob/{commit}{path}#{line}",
    sitename="Hestia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stephans3.github.io/Hestia.jl",
        assets=String[],
    ),
    pages=Any[
       "Hestia.jl" => "index.md",
       "Theory" => Any[
         "theory/geometry_material.md",
         "theory/boundary_conditions.md",
         "theory/actuators.md"
         ],
       "Tutorials" => Any[
         "examples/rod_1d.md",
         "examples/plate_2d.md"
        ], 
   	    "Components" => "components.md"
    ],
)


deploydocs(;
    repo="github.com/stephans3/Hestia.jl.git",
)
