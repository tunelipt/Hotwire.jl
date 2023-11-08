using Hotwire
using Documenter

DocMeta.setdocmeta!(Hotwire, :DocTestSetup, :(using Hotwire); recursive=true)

makedocs(;
    modules=[Hotwire],
    authors="Paulo José Saiz Jabardo",
    repo="https://github.com/tunelipt/Hotwire.jl/blob/{commit}{path}#{line}",
    sitename="Hotwire.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=:doctest     
)
