using QEDprobing
using Documenter

DocMeta.setdocmeta!(QEDprobing, :DocTestSetup, :(using QEDprobing); recursive=true)

makedocs(;
    modules=[QEDprobing],
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de",
    sitename="QEDprobing.jl",
    format=Documenter.HTML(;
        canonical="https://qedapplications.gitlab.io/QEDprobing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
