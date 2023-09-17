using SCPModelingToolkit
using Documenter

DocMeta.setdocmeta!(SCPModelingToolkit, :DocTestSetup, :(using SCPModelingToolkit); recursive=true)

makedocs(;
    modules=[SCPModelingToolkit],
    authors="Ben Chung <ckfinite@gmail.com> and contributors",
    repo="https://github.com/benchung/SCPModelingToolkit.jl/blob/{commit}{path}#{line}",
    sitename="SCPModelingToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
