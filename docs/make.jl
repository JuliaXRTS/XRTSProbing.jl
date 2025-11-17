using Pkg

# targeting the correct source code
# this assumes the make.jl script is located in QEDbase.jl/docs
project_path = Base.Filesystem.joinpath(Base.Filesystem.dirname(Base.source_path()), "..")

using XRTSProbing
using Documenter

DocMeta.setdocmeta!(XRTSProbing, :DocTestSetup, :(using XRTSProbing); recursive = true)

# some paths for links
readme_path = joinpath(project_path, "README.md")
index_path = joinpath(project_path, "docs/src/index.md")
license_path = "https://github.com/QEDjl-project/QEDbase.jl/blob/main/LICENSE"
coc_path = "https://github.com/QEDjl-project/QEDbase.jl/blob/main/CODE_OF_CONDUCT.md"

# Copy README.md from the project base folder and use it as the start page
open(readme_path, "r") do readme_in
    readme_string = read(readme_in, String)

    # replace relative links in the README.md
    readme_string = replace(readme_string, "[MIT](LICENSE)" => "[MIT]($(license_path))")
    readme_string = replace(readme_string, "(docs/src/90-contributing.md)" => "(90-contributing.md)")
    readme_string = replace(readme_string, "(docs/src/90-contributing.md)" => "(90-contributing.md)")
    open(index_path, "w") do readme_out
        write(readme_out, readme_string)
    end
end

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
        file != "index.md" && splitext(file)[2] == ".md"
]

try
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

finally
    # doing some garbage collection
    @info "GarbageCollection: remove generated landing page"
    rm(index_path)
end

deploydocs(; repo = "github.com/JuliaXRTS/XRTSProbing.jl.git", push_preview = true)
