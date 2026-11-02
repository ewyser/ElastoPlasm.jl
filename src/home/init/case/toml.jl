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
    save_case(file, case) -> file
    save_case(file, case, solver; folder=case_folder(case)) -> file

Write a case to a TOML file: its layout with its settings (EasyDraw's `save_layout`), SI throughout (φ in radians):

- the geometry: `[domain]` (`L`, `nel`, `ni`) and one `[[zone]]` per zone (its name and shape: `shape = "polygon"`
  with `vertices`, or `rings` for several outlines; `"block"` with `lo`/`hi`; `"ball"` with `centre`/`radius`);
- the settings: `[settings.zones.<name>]` (each zone's material constants or velocity `v`), `[settings.defaults]` (the
  material constants set), `[settings.sides]` (a boundary condition per side: `left`, `right`, `bottom`, `top`, and in
  3D `front`, `back`), `[settings.options.<section>]` (the solver overrides: numbers and booleans as such, other
  values as their `get_option()` label, `option_labels`; quantities by key), `[settings.run]` (`folder`, `te`, `tg`,
  `tep`, `workflow`, `lstps`, and `thermal` when true).

With a solver, the case as that solver solves it (what `setup_simulation` writes beside a simulation): the numerical
options are the solver's (those that differ from `get_default()`), the run settings its own and `folder`, so that
`Solver(load_case(file))` is the same solver.

A zone given as a Julia predicate (`Region`) cannot be written: an `ArgumentError`. To keep geometry and settings
apart, write them with `save_layout(file, case.layout)` and `save_settings(file, case.settings)`.
"""
function save_case(file::AbstractString, case::Case)
    return save_layout(file, case.layout, case.settings)
end
function save_case(file::AbstractString, case::Case, solver::Solver; folder::AbstractString=case_folder(case))
    s = case.settings
    c = solver_config(solver)
    numerics = config_changes(Base.structdiff(c, NamedTuple{(:run,)}), get_default())
    options = Dict{String,Any}(k => v for (k, v) ∈ s.options if Symbol(k) ∈ get_physics_sections())
    merge!(options, option_labels(config_dict(numerics)))
    run = merge(Dict{String,Any}(s.run), config_dict(c.run), Dict{String,Any}("folder" => folder))
    settings = EasyDraw.Settings(; zones=s.zones, defaults=s.defaults, sides=s.sides, options, run)
    return save_layout(file, case.layout, settings)
end

"""
    load_case(file) -> Case
    load_case(geometry, settings) -> Case

Read a case from its file (`save_case`): any layout file with settings whose values ElastoPlasm knows, so also one
drawn in EasyDraw's window with `get_catalog()`; or from a geometry file and a settings file, so that one geometry
serves several settings. The run folder is the settings' `[run]` folder, else the geometry file's own folder.

A case file of ElastoPlasm's older format (with `bcs`, `[material]`, `[solver]`, and φ in degrees) is refused, with an
`ArgumentError`: set the case up again from its example, or as `Case(; …)`, and save it again.
"""
function load_case(file::AbstractString)
    return load_case(file, file)
end
function load_case(geometry::AbstractString, settings::AbstractString)
    for file ∈ unique([geometry, settings])
        old = intersect(keys(TOML.parsefile(file)), ("bcs", "material", "solver"))
        isempty(old) || throw(ArgumentError("$(basename(file)) is a case file of the older format ($(join(sort(collect(old)), ", "))): " *
                                            "set the case up again from its example, or as Case(; …), and save it again"))
    end
    return Case(load_layout(geometry), load_settings(settings))
end
