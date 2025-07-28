using Printf, REPL.TerminalMenus

function read_version(path = joinpath(@__DIR__, "..", "Project.toml"))
    for line in eachline(path)
        if occursin("version", line)
            return match(r"\"(.*?)\"", line).captures[1]
        end
    end
    error("Version not found in Project.toml")
end

function write_version(new_version; path = joinpath(@__DIR__, "..", "Project.toml"))
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
println("Current version: v$current_version")

# Step 2: bump patch version
part = select_version_part()
new_version = bump_version(current_version, part)
println("Bumping $part version to: v$new_version")

# Step 3: write new version back to Project.toml
write_version(new_version)

# Step 4: git add Project.toml
run(`git add ../Project.toml`)

# Step 5: git commit with new version
run(`git commit -m "Release v$new_version"`)

# Step 6: git tag
run(`git tag v$new_version`)

# Step 7: push commit and tag
run(`git push origin main`)
run(`git push origin v$new_version`)