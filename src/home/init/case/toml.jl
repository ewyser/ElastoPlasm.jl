export load_case, save_case

# ---- solver options on disk: values in memory, labels or plain numbers in files

# the `(prompt, values)` entry of `get_option()` at `path` (e.g. ["basis", "which"]), or nothing
function option_entry(path)
    o = get_option()
    for k ∈ path
        (o isa NamedTuple && haskey(o, Symbol(k))) || return nothing
        o = getfield(o, Symbol(k))
    end
    return o isa Tuple && length(o) == 2 && o[1] isa AbstractString ? o : nothing
end
option_label(v) = string(v)
option_label(v::NamedTuple) = haskey(v, :mpts) ? v.mpts.name : v.mesh.name        # a plot.what variable: its name
# the `get_default()` value at `path`
default_at(path) = foldl((d, k) -> getfield(d, Symbol(k)), path; init=get_default())
# a numeric option takes any value, not only the listed ones
free_numeric(d) = (d isa Real && !(d isa Bool)) || d isa AbstractVector{<:Real}
# a multi-selection: an option whose value is a list of the listed values (e.g. `plot.what`)
multi_select(path) = (d = default_at(path); d isa AbstractVector && !(d isa AbstractVector{<:Real}))
# solver options: values in memory; on disk booleans, numbers and number lists as such, the rest as their
# `get_option()` labels (so e.g. `dtype` types round-trip), plot variables by name
option_save(v) = v isa Real || v isa AbstractVector{<:Real} ? v : v isa AbstractVector ? option_label.(v) : option_label(v)
option_labels(d::AbstractDict) = Dict{String,Any}(k => v isa AbstractDict ? option_labels(v) : option_save(v) for (k, v) ∈ d)
option_values(d::AbstractDict, path=String[]) =
    Dict{String,Any}(k => v isa AbstractDict ? option_values(v, [path; k]) : option_value([path; k], v) for (k, v) ∈ d)
# a number list typed as text: "2.5, 2.5, 2.5" or "[2.5, 2.5, 2.5]"
parse_numbers(s::AbstractString) = (t = split(strip(s, ['[', ']', ' ']), r"[,\s]+"; keepempty=false);
                                    v = tryparse.(Float64, t); isempty(v) || any(isnothing, v) ? nothing : Float64.(v))
"""
    option_value(path, x) -> value

The solver option at `path` (e.g. `["basis", "which"]`) for `x`, a value or its label: one of the
`get_option()` values, any number (list) for a numeric option (converted to the default's type), and for
a multi-selection (e.g. `plot.what`) a list of the listed values or their names.
"""
function option_value(path, x)
    e = option_entry(path)
    isnothing(e) && throw(ArgumentError("unknown solver option $(join(path, "."))"))
    vals, name = e[2], join(path, ".")
    if multi_select(path)
        x isa AbstractVector || throw(ArgumentError("solver option $name is a list of $(option_label.(vals))"))
        return [option_value_in(vals, s, name) for s ∈ x]
    end
    i = findfirst(v -> option_label(v) == option_label(x), vals)
    isnothing(i) || return vals[i]
    d = default_at(path)
    free_numeric(d) || throw(ArgumentError("solver option $name = $x is not one of $(option_label.(vals))"))
    y = x isa AbstractString ? (d isa AbstractVector ? parse_numbers(x) : tryparse(Float64, x)) : x
    if d isa AbstractVector
        y isa AbstractVector{<:Real} || throw(ArgumentError("solver option $name = $x is not a list of numbers"))
        return Float64.(y)
    end
    y isa Real && !(y isa Bool) || throw(ArgumentError("solver option $name = $x is not a number"))
    d isa Integer && !isinteger(y) && throw(ArgumentError("solver option $name = $x is not an integer"))
    return d isa Integer ? typeof(d)(y) : Float64(y)
end
option_value_in(vals, s, name) = (i = findfirst(v -> option_label(v) == string(s), vals);
    isnothing(i) ? throw(ArgumentError("solver option $name: $s is not one of $(option_label.(vals))")) : vals[i])

# a NamedTuple (or Dict) of option sections as a Dict with string keys, and back
config_dict(x) = x isa NamedTuple || x isa AbstractDict ? Dict{String,Any}(string(k) => config_dict(v) for (k, v) ∈ pairs(x)) : x

# ---- the case on disk

const SHAPES = ("polygon", "block", "ball")
deg(k, v) = k == :phi ? round(rad2deg(v); sigdigits=12) : v   # angles: degrees in files (12 digits: 30, not 29.999999999999996)
rad(k, v) = k == :phi ? deg2rad(v) : v

function shape_dict(s::AbstractShape)
    s isa Polygon && return Dict{String,Any}("shape" => "polygon", "vertices" => s.vertices)
    s isa Block     && return Dict{String,Any}("shape" => "block", "lo" => s.lo, "hi" => s.hi)
    s isa Ball    && return Dict{String,Any}("shape" => "ball", "centre" => s.centre, "radius" => s.radius)
    throw(ArgumentError("a $(nameof(typeof(s))) zone (a Julia predicate) cannot be written to TOML"))
end
function parse_shape(z::AbstractDict, name)
    kind = get(z, "shape", "")
    kind == "polygon" && return Polygon(z["vertices"])
    kind == "block"     && return Block(z["lo"], z["hi"])
    kind == "ball"    && return Ball(z["centre"], z["radius"])
    throw(ArgumentError("zone \"$name\": shape = \"$kind\" is not one of $(SHAPES)"))
end

"""
    case_dict(case) -> Dict

The TOML form of a case (`save_case`): `[domain]`, the material constants differing from
`MATERIAL`, `bcs`, `thermal`, `[solver.…]` overrides (`option_save`), `[run]` and one `[[zone]]` per
zone; angles in degrees.
"""
function case_dict(s::Case)
    dom = Dict{String,Any}("L" => s.L, "nel" => s.nel, "ni" => s.ni)
    isnothing(s.width) || (dom["width"] = s.width; dom["nely"] = s.nely)
    zones = map(s.zones) do z
        d = shape_dict(z.shape)
        d["name"] = z.name
        for (k, v) ∈ z.props
            d[string(k)] = deg(k, v)
        end
        d
    end
    return Dict{String,Any}("domain" => dom,
        "material" => Dict{String,Any}(string(k) => deg(k, v) for (k, v) ∈ pairs(s.material) if v != getfield(MATERIAL, k)),
        "bcs" => [string.(s.bcs[d, :]) for d ∈ axes(s.bcs, 1)], "thermal" => s.thermal,
        "solver" => option_labels(config_dict(s.solver)), "run" => config_dict(s.run), "zone" => zones)
end

"""
    parse_case(d) -> Case

A problem from its TOML form (`case_dict`, `load_case`).
"""
function parse_case(d::AbstractDict)
    haskey(d, "domain") || throw(ArgumentError("a case file needs a [domain] table (L, nel, ni)"))
    dom   = d["domain"]
    zones = map(enumerate(get(d, "zone", []))) do (i, z)
        name = String(get(z, "name", "zone $i"))
        known = ("name", "shape", "vertices", "lo", "hi", "centre", "radius")
        Zone(name, parse_shape(z, name), Dict{Symbol,Any}(Symbol(k) => rad(Symbol(k), v) for (k, v) ∈ z if k ∉ known))
    end
    mat = Dict{Symbol,Any}(Symbol(k) => rad(Symbol(k), v) for (k, v) ∈ get(d, "material", Dict()))
    bcs = haskey(d, "bcs") ? permutedims(reduce(hcat, [Symbol.(side) for side ∈ d["bcs"]])) : nothing
    return Case(; L=dom["L"], nel=dom["nel"], ni=get(dom, "ni", 2), width=get(dom, "width", nothing),
                  nely=get(dom, "nely", nothing), zones, material=mat, bcs, thermal=get(d, "thermal", false),
                  solver=option_values(get(d, "solver", Dict())), run=get(d, "run", Dict()))
end

"""
    save_case(file, case) -> file
    load_case(file) -> Case

Write a case to a TOML file, and read it back: `[domain]` (`L`, `nel`, `ni`, optionally `width` and
`nely`), `[material]` (overrides of `MATERIAL`, `phi` in degrees), `bcs` (`[[lower, upper], …]` per
axis), `thermal`, `[solver.<section>]` overrides (values as `get_option()` lists them, any number for
a numeric option, variable names for `plot.what`), `[run]` (`folder`, `te`, `tg`, `tep`, `workflow`,
`lstps`) and `[[zone]]` tables (`name`, `shape = "polygon"` with `vertices = [[x, z], …]`, `"block"` with
`lo`/`hi`, or `"ball"` with `centre`/`radius`, and any material constant or `v`). A zone given as a
Julia predicate (`Region`) cannot be written.
"""
function save_case(file::AbstractString, s::Case)
    d = case_dict(s)
    mkpath(dirname(abspath(file)))
    open(io -> TOML.print(io, d; sorted=true), file, "w")
    return file
end
load_case(file::AbstractString) = parse_case(TOML.parsefile(file))
