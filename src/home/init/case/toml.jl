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
# the `get_default()` value at `path`
default_at(path) = foldl((d, k) -> getfield(d, Symbol(k)), path; init=get_default())
# a numeric option takes any value, not only the listed ones
free_numeric(d) = (d isa Real && !(d isa Bool)) || d isa AbstractVector{<:Real}
# a multi-selection: an option whose value is a list of the listed values (`plot.what`, the `output.what` lists)
multi_select(path) = (d = default_at(path); d isa AbstractVector && !(d isa AbstractVector{<:Real}))
# solver options: values in memory; on disk booleans, numbers and number lists as such, the rest as their
# `get_option()` labels (so e.g. `dtype` types round-trip), quantities by key
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
a multi-selection (`plot.what`, the `output.what` lists) a list of the listed values or their names.
"""
function option_value(path, x)
    e = option_entry(path)
    isnothing(e) && throw(ArgumentError("unknown solver option $(join(path, "."))"))
    vals, name = e[2], join(path, ".")
    if multi_select(path)
        x isa AbstractVector || throw(ArgumentError("solver option $name is a list of $(option_label.(vals))"))
        return [option_value_in(vals, s, path) for s ∈ x]
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
    y isa Real && !(y isa Bool) && isfinite(y) || throw(ArgumentError("solver option $name = $x is not a (finite) number"))
    d isa Integer && !isinteger(y) && throw(ArgumentError("solver option $name = $x is not an integer"))
    return d isa Integer ? typeof(d)(y) : Float64(y)
end
# a value of a multi-selection, by label or by its quantity's readable name (older cases wrote `plot.what` so)
function option_value_in(vals, s, path)
    i = something(findfirst(v -> option_label(v) == string(s), vals), findfirst(v -> shown_name(path, v) == string(s), vals), 0)
    i == 0 && throw(ArgumentError("solver option $(join(path, ".")): $s is not one of $(option_label.(vals))"))
    return vals[i]
end

# a NamedTuple (or Dict) of option sections as a Dict with string keys, and back
config_dict(x) = x isa NamedTuple || x isa AbstractDict ? Dict{String,Any}(string(k) => config_dict(v) for (k, v) ∈ pairs(x)) : x

# ---- the case on disk: its layout's and settings' file

"""
    save_case(file, case, solver; folder=case_folder(case)) -> file

Write a case as `solver` solves it (what `setup` writes beside a simulation): EasyDraw's `save_case`
(the layout with its settings, SI throughout, φ in radians), with the settings' numerical options replaced by the
solver's (those that differ from `get_default()`, as their `get_option()` labels, `option_labels`) and its run settings
with `folder`, so that `Solver(load_case(file))` is the same solver. The physics options (`material`, `grf`) are the
case's.

A case file holds, under `[settings]`: `[settings.zones.<name>]` (each zone's material constants or velocity `v`),
`[settings.defaults]` (the material constants set), `[settings.sides]` (a boundary condition per side, roller where
not set), `[settings.options.<section>]`, `[settings.run]` (`folder`, `te`, `tg`, `tep`, `workflow`, `lstps`, and
`thermal` when true). EasyDraw's `load_case` reads it back, and refuses a file of ElastoPlasm's older format (`bcs`,
`[material]`, `[solver]`, φ in degrees): set such a case up again from its example.

# Throws
- `ArgumentError`: when a zone is a `Region` predicate, which a file cannot hold.
"""
function EasyDraw.save_case(file::AbstractString, case::Case, solver::Solver; folder::AbstractString=case_folder(case))
    s = case.settings
    c = solver_config(solver)
    numerics = config_changes(Base.structdiff(c, NamedTuple{(:run,)}), get_default())
    options = Dict{String,Any}(k => v for (k, v) ∈ s.options if Symbol(k) ∈ get_physics_sections())
    merge!(options, option_labels(config_dict(numerics)))
    run = merge(Dict{String,Any}(s.run), config_dict(c.run), Dict{String,Any}("folder" => folder))
    settings = EasyDraw.Settings(; zones=s.zones, defaults=s.defaults, sides=s.sides, options, run)
    return save_case(file, Case(case.layout, settings))
end
