export get_catalog, write_catalog                            # EasyDraw's write_catalog: get_catalog() as a TOML file

# the options pages' order in the window (a section not listed comes last)
const PAGE_ORDER = ["basis", "stab", "material", "grf", "nonloc", "plot", "output", "dtype", "backend"]

"""
    get_catalog() -> Dict

What the drawing window offers for a case (`draw_case`), as plain data shaped like a catalog's TOML (EasyDraw reads
it: `EasyDraw.Catalog`, `write_catalog`), built from the tables that define it: the zone constants (c₀, cᵣ, φ, T, E,
ν, ρ₀; their defaults `MATERIAL`, φ only under Drucker–Prager), the boundary kinds (`BOUNDARY_KINDS`, roller by
default), every solver option of `get_option()` with its `get_default()` value (but `dtype.T0`/`dtype.bits`, set from
`dtype.precision` by `get_solver`, and `solution`, set by the workflow), quantities by their readable names
(`shown_name`), and the run settings (`WORKFLOWS`; `te`, `tg`, `tep` for the explicit workflows, `lstps` for the
quasi-static ones; `RUN`'s defaults). One `Dict` per field: `key`, `kind` (`"menu"`, `"numeric"`, `"integer"`,
`"list"`), `label`, `values`, `default`, `scale` (shown = SI × scale), `min`, `show_if`, `names`.
"""
function get_catalog()
    zone = [catalog_field("c0"; label="c₀ [kPa]", scale=1e-3, default=MATERIAL.c0),
            catalog_field("cr"; label="cᵣ [kPa]", scale=1e-3, default=MATERIAL.cr),
            catalog_field("phi"; label="φ [°]", scale=180/π, default=MATERIAL.phi, show_if=Dict("material.plastic" => ["DP"])),
            catalog_field("T"; label="T [K]", default=MATERIAL.T), catalog_field("E"; label="E [MPa]", scale=1e-6, default=MATERIAL.E),
            catalog_field("nu"; label="ν [–]", default=MATERIAL.nu), catalog_field("rho0"; label="ρ₀ [kg/m³]", default=MATERIAL.rho0)]
    options = Dict{String,Any}[]
    function walk(o, d, sec, path)
        if o isa Tuple && length(o) == 2 && o[1] isa AbstractString
            skip = (sec == "dtype" && path ∈ (["T0"], ["bits"])) || sec == "solution" || isempty(o[2])
            skip || push!(options, option_field(sec, path, o[2], d))
        elseif o isa NamedTuple
            foreach(k -> walk(getfield(o, k), getfield(d, k), sec, [path; string(k)]), keys(o))
        end
    end
    opt, def = get_option(), get_default()
    foreach(k -> walk(getfield(opt, k), getfield(def, k), string(k), String[]), keys(opt))
    page(f) = something(findfirst(==(first(split(f["key"], "."))), PAGE_ORDER), length(PAGE_ORDER) + 1)
    explicit, quasi = [w for w ∈ WORKFLOWS if !quasistatic(w)], [w for w ∈ WORKFLOWS if quasistatic(w)]
    run = [catalog_field("workflow"; kind="menu", label="workflow", values=collect(WORKFLOWS), default=RUN.workflow),
           catalog_field("te"; label="elastic stage te [s]", default=RUN.te, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("tg"; label="gravity ramp tg [s]", default=RUN.tg, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("tep"; label="plastic stage tep [s]", default=RUN.tep, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("lstps"; kind="integer", label="load steps", default=RUN.lstps, min=1, show_if=Dict("workflow" => quasi))]
    return Dict{String,Any}("zone" => zone, "options" => sort(options; by=page), "run" => run,
                            "boundary" => Dict{String,Any}("kinds" => string.(collect(BOUNDARY_KINDS)),
                                                           "colors" => ["#cccccc", "#4682b4", "#b22222"], "default" => "roller"))
end
# one field of the catalog: its key and the entries given (an entry left `nothing` is left out)
function catalog_field(key; kwargs...)
    return Dict{String,Any}("key" => key, (string(k) => v for (k, v) ∈ kwargs if !isnothing(v))...)
end
# the option at `path` below section `sec`: a list (`multi_select`, labelled by its group, "point · solid", its values by
# their readable names), a number (any value, `free_numeric`), or a drop-down
function option_field(sec, path, vals, d)
    key = join([sec; path], ".")
    if multi_select([sec; path])
        return catalog_field(key; kind="list", label=length(path) > 1 ? join(path[2:end], " · ") : join(path, "."), values=vals,
                             default=d, names=[shown_name([sec; path], v) for v ∈ vals])
    end
    kind = free_numeric(d) ? (d isa Integer ? "integer" : "numeric") : "menu"
    return catalog_field(key; kind, label=join(path, "."), values=vals, default=d)
end
