using XRTSProbing
using Documenter

DocMeta.setdocmeta!(XRTSProbing, :DocTestSetup, :(using XRTSProbing); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
        file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [XRTSProbing],
    authors = "Uwe Hernandez Acosta <u.hernandez@hzdr.de>",
    repo = "https://github.com/JuliaXRTS/XRTSProbing.jl/blob/{commit}{path}#{line}",
    sitename = "XRTSProbing.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://github.com/JuliaXRTS/XRTSProbing.jl/blob/{commit}{path}#{line}",
        canonical = "https://JuliaXRTS.github.io/XRTSProbing.jl",
        assets = String[],
    ),
    pages = ["index.md"; numbered_pages],
)

#deploydocs(; repo = "github.com/QEDjl-project/QEDprobing.jl")
