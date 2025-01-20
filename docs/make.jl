using QEDprobing
using Documenter

DocMeta.setdocmeta!(QEDprobing, :DocTestSetup, :(using QEDprobing); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [QEDprobing],
    authors = "Uwe Hernandez Acosta <u.hernandez@hzdr.de>",
    repo = "https://codebase.helmholtz.cloud/qedjl-applications/QEDprobing.jl/blob/{commit}{path}#{line}",
    #repo = "https://github.com/QEDjl-project/QEDprobing.jl/blob/{commit}{path}#{line}",
    sitename = "QEDprobing.jl",
    #format = Documenter.HTML(; canonical = "https://QEDjl-project.github.io/QEDprobing.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://codebase.helmholtz.cloud/qedjl-applications/QEDprobing.jl/blob/{commit}{path}#{line}",
        canonical = "https://qedjl-applications.pages.hzdr.de/qedprobing.jl",
        assets = String[],
    ),
    pages = ["index.md"; numbered_pages],
)

#deploydocs(; repo = "github.com/QEDjl-project/QEDprobing.jl")
