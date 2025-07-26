using Documenter,ElastoPlasm

DocMeta.setdocmeta!(ElastoPlasm, :DocTestSetup, :(using ElastoPlasm); recursive=true)

format = Documenter.HTML()
manual = [
    "Home"            => "index.md",
    "Getting Started" => "getting_started.md",
    #"Functions" => mdGenerate(),
]

repo = "https://github.com/ewyser/ElastoPlasm.jl.git"

@info "Making documentation..."
makedocs(;
    modules=[ElastoPlasm],
    authors="madmax",
    sitename="ϵlastσPlasm.jl 👻",
    format=Documenter.HTML(;
        repolink = repo,
        canonical = "https://ewyser.github.io/ElastoPlasm.jl/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Function" => "function.md",
    ],
    checkdocs= :none,
)

@info "Deploying documentation..."

allowed_branches = ["refs/heads/main", "refs/heads/misc"]
current_ref = get(ENV, "GITHUB_REF", "")

if current_ref ∈ allowed_branches
    deploydocs(;repo = repo,
                devbranch = "main",
                branch = "gh-pages",
                versions = ["stable" => "v^", "dev" => "dev"],
                forcepush = true,
                push_preview = true,
            )
else
    @info "Current branch $current_ref is not in allowed branches $allowed_branches, skipping deploydocs."
end
#=
@info "Deploying documentation..."
if get(ENV, "GITHUB_EVENT_NAME", "") == "pull_request"

else
    withenv("GITHUB_REPOSITORY" => repo) do
        deploydocs(;repo = repo,
                    devbranch = "main",
                    branch = "gh-pages",
                    versions = nothing,#["stable" => "v^", "dev" => "dev", "v#.#.#"],
                    forcepush = true,
                    push_preview = true,
                )
    end
end
=#