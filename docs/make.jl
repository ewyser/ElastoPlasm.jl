
using Documenter, ElastoPlasm

DocMeta.setdocmeta!(ElastoPlasm, :DocTestSetup, :(using ElastoPlasm); recursive=true)

repo = "https://github.com/ewyser/ElastoPlasm.jl.git"

# Helper function to generate an @autodocs block string
function write_autodoc_page(filename, pages, modulename)
    check = joinpath(@__DIR__, "src","function")
    if !isdir(check)
        mkpath(check)    
    end
    path = joinpath(check, filename)
    
    open(path, "w") do io
        println(io, "# $modulename Reference\n")
        println(io, "```@meta")
        println(io, "CollapsedDocStrings = true")
        println(io, "```")
        println(io)
        println(io, "```@autodocs")
        println(io, "Modules = [$modulename]")
        println(io, "Pages   = ", repr(pages))
        println(io, "```")
    end
end
# Recursively collect .jl filenames from a superInc nested dict
function collect_jl_files(d::Dict, path::String = "")
    files = String[]
    for (k, v) in d
        if isa(v, Dict)
            append!(files, collect_jl_files(v, joinpath(path, k)))
        elseif endswith(k, ".jl")
            push!(files, isempty(path) ? k : joinpath(path, k))
        end
    end
    return files
end

lib_pages(key) = collect_jl_files(ElastoPlasm.self.sys.lib[key]["files"])

write_autodoc_page("api.md"    , lib_pages("home/api")   , :ElastoPlasm)
write_autodoc_page("program.md", lib_pages("home/core")  , :ElastoPlasm)
write_autodoc_page("script.md" , lib_pages("home/script"), :ElastoPlasm)

# Call makedocs
@info "Making documentation..."
makedocs(;
    modules = [ElastoPlasm],
    authors = "madmax",
    sitename = "ϵlastσPlasm.jl 👻",
    format = Documenter.HTML(;
        repolink = repo,
        canonical = "https://ewyser.github.io/ElastoPlasm.jl/",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Function Reference" => [
            "Public API" => "function/api.md",
            "Internals"  => "function/program.md",
            "Example"    => "function/script.md",
        ]
    ],
    checkdocs = :none,
)

# Deploy documentation only if inside GitHub Actions
if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    @info "Deploying documentation..."
    deploydocs(; 
        repo = repo,
        devbranch = "main",
        branch = "gh-pages",
        versions = ["stable" => "v^", "dev" => "dev"],
        forcepush = true,
        push_preview = true,
    )
else
    @info "Not running inside CI, skipping deploydocs."
end