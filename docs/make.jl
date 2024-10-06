using EigenMicrostates
using Documenter

DocMeta.setdocmeta!(EigenMicrostates, :DocTestSetup, :(using EigenMicrostates); recursive=true)

makedocs(;
    modules=[EigenMicrostates],
    authors="Zanchenling Wang <43379393+wangzcl@users.noreply.github.com> and contributors",
    sitename="EigenMicrostates.jl",
    format=Documenter.HTML(;
        canonical="https://wangzcl.github.io/EigenMicrostates.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wangzcl/EigenMicrostates.jl",
    devbranch="main",
)
