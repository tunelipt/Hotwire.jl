using Lixo
using Documenter

DocMeta.setdocmeta!(Lixo, :DocTestSetup, :(using Lixo); recursive=true)

makedocs(;
    modules=[Lixo],
    authors="Paulo José Saiz Jabardo",
    repo="https://github.com/pjsjipt/Lixo.jl/blob/{commit}{path}#{line}",
    sitename="Lixo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
