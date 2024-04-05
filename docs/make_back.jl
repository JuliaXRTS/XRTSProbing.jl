import Pkg
using QEDprobing
using Documenter
using DocumenterMermaid

# targeting the correct source code
# this asumes the make.jl script is located in QEDbase.jl/docs

project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")
Pkg.develop(; path=project_path)
DocMeta.setdocmeta!(QEDprobing, :DocTestSetup, :(using QEDprobing); recursive=true)

makedocs(;
    modules=[QEDprobing],
         checkdocs=:exports,
    authors="Uwe Hernandez Acosta <u.hernandez@hzdr.de",
    sitename="QEDprobing.jl",
    repo=Documenter.Remotes.GitLab("codebase.helmholtz.cloud","qedjl-applications","QEDprobing.jl"),
    format=Documenter.HTML(;
                           prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://qedapplications.gitlab.io/QEDprobing.jl",
        #edit_link="main",
        assets=String[],
    ),
         pages=[
         "Home" => "index.md",
         "Sampling Algorithms" => [
           "sampling_algorithms/overview.md",
           "sampling_algorithms/event_generation.md",
           "sampling_algorithms/compound_distribution.md",
         ],
         ],
         )
