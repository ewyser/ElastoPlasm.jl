using Documenter,ElastoPlasm

DocMeta.setdocmeta!(ElastoPlasm, :DocTestSetup, :(using ElastoPlasm); recursive=true)

format = Documenter.HTML()
manual = [
    "Home"            => "index.md",
    "Getting Started" => "getting_started.md",
    #"Functions" => mdGenerate(),
]
@info "Making documentation..."
makedocs(;
    modules=[ElastoPlasm],
    authors="madmax",
    sitename="ϵlastσPlasm.jl 👻",
    format=Documenter.HTML(;
        repolink="github.com/ewyser/ElastoPlasm.jl",
        canonical="https://ewyser.github.io/ElastoPlasm.jl/",
        edit_link="misc",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs= :none,
)

@info "Deploying documentation..."

allowed_branches = ["refs/heads/main", "refs/heads/misc"]

current_ref = get(ENV, "GITHUB_REF", "")

if !(current_ref ∈ allowed_branches)
    @info "Current branch $current_ref is not in allowed branches $allowed_branches, skipping deploydocs."
    return
end

deploydocs(
    repo = "https://github.com/ewyser/ElastoPlasm.jl.git",
    branch = "gh-pages",
    target = ""#"misc",
)
