export get_catalog, write_catalog                            # EasyDraw's write_catalog: get_catalog() as a TOML file

"""
    get_catalog() -> Dict

What the drawing window offers for a case (`draw_case`), as plain data shaped like a catalog's TOML (EasyDraw reads it:
`EasyDraw.Catalog`, `write_catalog`), built from the tables that define it: the zone constants (c₀, cᵣ, φ, T, E, ν, ρ₀;
their defaults `get_material()`, φ only under Drucker–Prager), the boundary kinds (`get_boundary_kinds()`, roller by
default), every solver option of `get_option()` with its `get_default()` value (but `dtype.T0`/`dtype.bits`, set from
`dtype.precision` by `solver_instr`, and `solution`, set by the workflow), quantities by their readable names
(`shown_name`), and the run settings (`get_workflows()`; `te`, `tg`, `tep` for the explicit workflows, `lstps` for the
quasi-static ones; `get_run()`'s defaults). One `Dict` per field: `key`, `kind` (`"menu"`, `"numeric"`, `"integer"`,
`"list"`), `label`, `values`, `default`, `scale` (shown = SI × scale), `min`, `show_if`, `names`.
"""
function get_catalog()
    m, r, w = get_material(), get_run(), get_workflows()
    zone = [catalog_field("c0"; label="c₀ [kPa]", scale=1e-3, default=m.c0),
            catalog_field("cr"; label="cᵣ [kPa]", scale=1e-3, default=m.cr),
            catalog_field("phi"; label="φ [°]", scale=180/π, default=m.phi, show_if=Dict("material.plastic" => ["DP"])),
            catalog_field("T"; label="T [K]", default=m.T), catalog_field("E"; label="E [MPa]", scale=1e-6, default=m.E),
            catalog_field("nu"; label="ν [–]", default=m.nu), catalog_field("rho0"; label="ρ₀ [kg/m³]", default=m.rho0)]
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
    pages = ["basis", "stab", "material", "grf", "nonloc", "plot", "output", "dtype", "backend"]   # the window's order
    page(f) = something(findfirst(==(first(split(f["key"], "."))), pages), length(pages) + 1)          # others last
    explicit, quasi = [x for x ∈ w if !quasistatic(x)], [x for x ∈ w if quasistatic(x)]
    run = [catalog_field("workflow"; kind="menu", label="workflow", values=collect(w), default=r.workflow),
           catalog_field("te"; label="elastic stage te [s]", default=r.te, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("tg"; label="gravity ramp tg [s]", default=r.tg, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("tep"; label="plastic stage tep [s]", default=r.tep, min=0, show_if=Dict("workflow" => explicit)),
           catalog_field("lstps"; kind="integer", label="load steps", default=r.lstps, min=1, show_if=Dict("workflow" => quasi))]
    return Dict{String,Any}("zone" => zone, "options" => sort(options; by=page), "run" => run,
                            "boundary" => Dict{String,Any}("kinds" => string.(collect(get_boundary_kinds())),
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
