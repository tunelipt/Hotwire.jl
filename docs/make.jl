using Hotwire
using Documenter
import Polynomials

DocMeta.setdocmeta!(Hotwire, :DocTestSetup, :(using Hotwire); recursive=true)

makedocs(;
    modules=[Hotwire],
    authors="Paulo José Saiz Jabardo",
#    repo="https://github.com/tunelipt/Hotwire.jl/blob/{commit}{path}#{line}",
    sitename="Hotwire.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/tunelipt/Hotwire.jl/blob/{commit}{path}#{line}",
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=:doctest     
)
