using RobustPade
using Documenter

DocMeta.setdocmeta!(RobustPade, :DocTestSetup, :(using RobustPade); recursive=true)

makedocs(;
    modules=[RobustPade],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/RobustPade.jl/blob/{commit}{path}#{line}",
    sitename="RobustPade.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/RobustPade.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/RobustPade.jl",
    devbranch="main",
)
