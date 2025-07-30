using Printf, REPL.TerminalMenus

elastoplasm_project = joinpath(@__DIR__, "..", "Project.toml")

function read_version(path = elastoplasm_project)
    for line in eachline(path)
        if occursin("version", line)
            return match(r"\"(.*?)\"", line).captures[1]
        end
    end
    error("Version not found in Project.toml")
end

function write_version(new_version; path = elastoplasm_project)
    content = read(path, String)
    new_content = replace(content, r"version\s*=\s*\"[0-9]+\.[0-9]+\.[0-9]+\"" => "version = \"$new_version\"")
    open(path, "w") do io
        write(io, new_content)
    end
end

function bump_version(version::AbstractString, part::Symbol)
    major, minor, patch = parse.(Int, split(version, "."))
    if part == :major
        major += 1
        minor = 0
        patch = 0
    elseif part == :minor
        minor += 1
        patch = 0
    elseif part == :patch
        patch += 1
    else
        error("Unknown part: $part")
    end
    return @sprintf("%d.%d.%d", major, minor, patch)
end

function select_version_part(; prompt = "Select version part to bump:")
    options  = ["major", "minor", "patch"]
    selected = request(prompt, RadioMenu(options))
    return Symbol(options[selected])
end

# Step 1: read current version
current_version = read_version()

# Step 2: bump {major|minor|patch} version
@info "Bumping version and Releasing"
part = select_version_part()
new_version = bump_version(current_version, part)

# Step 3: repare commit message, confirm commit and tag
commit_msg = "Bump $current_version → $new_version and Release v$new_version"
options    = ["yes", "no"]
selected   = request("$commit_msg ?", RadioMenu(options))
if options[selected] == "yes"
    println("Bumping Project.toml, commit and tag release...")

    # Step 3a: bump Project.toml
    write_version(new_version)

    # Step 3b: git add Project.toml
    run(`git add $elastoplasm_project`)

    # Step 3c: git commit with new version
    run(`git commit -m "Bump $current_version and release v$new_version"`) →

    # Step 3d: git tag
    run(`git tag v$new_version`)

    # Step 3e: push commit and tag
    run(`git push origin main`)
    run(`git push origin v$new_version`)
    #=    =#
else
    println("Aborting release process.")
    return
end
