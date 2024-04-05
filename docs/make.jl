using QEDprobing
using Documenter

DocMeta.setdocmeta!(QEDprobing, :DocTestSetup, :(using QEDprobing); recursive=true)

makedocs(;
    modules=[QEDprobing],
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de",
    sitename="QEDprobing.jl",
    repo=Documenter.Remotes.URL(
        "https://codebase.helmholtz.cloud/qedjl-applications/QEDprobing.jl/lblob/{commit}{path}#{line}",
    ),
    format=Documenter.HTML(;
        repolink="https://codebase.helmholtz.cloud/qedjl-applications/QEDprobing.jl/lblob/{commit}{path}#{line}",
        canonical="https://qedapplications.gitlab.io/QEDprobing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)
