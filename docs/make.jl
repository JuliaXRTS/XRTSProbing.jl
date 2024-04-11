using QEDprobing
using Documenter

DocMeta.setdocmeta!(QEDprobing, :DocTestSetup, :(using QEDprobing); recursive=true)

makedocs(;
    modules=[QEDprobing],
    checkdocs=:exports,
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de",
    sitename="QEDprobing.jl",
    repo=Documenter.Remotes.GitLab(
        "codebase.helmholtz.cloud", "qedjl-applications", "QEDprobing.jl"
    ),
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://qedapplications.gitlab.io/QEDprobing.jl",
        #edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)
