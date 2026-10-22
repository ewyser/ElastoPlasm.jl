export quantities

# a vector of D-component vectors as a D × n matrix, one column each (what HDF5 stores)
function columns(v)
    return reduce(hcat, v)
end

"""
    qty(data, name, unit="-"; label=nothing, display=(unit, 1.0), cb=:viridis, cblim=nothing) -> NamedTuple

One entry of `quantities()`: returns its arguments as a NamedTuple, the defaults filled in.
`data` gives a frame of the quantity in SI, `name` is readable, `unit` the SI unit (`"-"`: dimensionless).
A `label` makes the quantity plottable, shown in `display = (unit, factor)` (by default its SI unit) with
colour map `cb` and limits `cblim` (in the displayed unit; `nothing`: the data's range).

```julia
qty(m -> m.s.ρ, "density", "kg/m^3")
# (; data = m -> m.s.ρ, name = "density", unit = "kg/m^3", label = nothing,
#    display = ("kg/m^3", 1.0), cb = :viridis, cblim = nothing)
```
"""
function qty(data, name, unit="-"; label=nothing, display=(unit, 1.0), cb=:viridis, cblim=nothing)
    return (; data, name, unit, label, display, cb, cblim)
end

"""
    quantities() -> NamedTuple

Every quantity a run can plot or record, shaped like `output.what`: object (`point`, `mesh`) → phase
(`base` for the point fields outside any phase, `solid`, `fluid`, `thermal`) → key → entry
`(; data, name, unit, label, display, cb, cblim)`. `data` gives one frame in SI (`unit`): a value per point
or node, or components × count. A point entry with a `label` can be plotted (`plot.what`), shown in
`display = (unit, factor)` with colour map `cb` and limits `cblim` (`nothing`: the data's range).
`output.what` and `plot.what` select keys, `get_option` lists them, an HDF5 output stores each field's
`name`, `unit` and `label` as attributes. To add a quantity, add an entry here.
"""
function quantities()
    return (;
        point = (;
            base = (;
                x     = qty(m -> columns(m.x), "position", "m"),
                J     = qty(m -> m.J, "deformation determinant"; label=L"J", cblim=(0.5, 1.5)),
                n     = qty(m -> m.n, "porosity"; label=L"n", cb=:vik, cblim=(-1.0, 1.0)),
                Omega = qty(m -> m.Ω, "volume", "m^D"),
                ell   = qty(m -> columns(m.ℓ), "domain", "m"; label=L"\ell"),   # plotted as the points' GIMP boxes
            ),
            solid = (;
                u     = qty(m -> columns(m.s.u), "displacement", "m"),
                v     = qty(m -> columns(m.s.v), "velocity", "m/s"),
                p     = qty(m -> getfield.(m.s.σᵢ, :p), "pressure", "Pa"; label=L"p", display=("kPa", 1e-3), cblim=(0.0, 150)),
                dev   = qty(m -> columns([σ.dev for σ ∈ m.s.σᵢ]), "deviatoric stress", "Pa"),
                epII  = qty(m -> getindex.(m.s.ϵpII, 1), "plastic strain";
                            label=L"\epsilon_{\mathrm{II}}^{\mathrm{acc}}", cblim=(0.0, 1.5)),
                epV   = qty(m -> m.s.ϵpV, "volumetric plastic strain"),
                rho   = qty(m -> m.s.ρ, "density", "kg/m^3"),
                q     = qty(m -> [sqrt(3*get_J2(σ)) for σ ∈ m.s.σᵢ], "von Mises stress", "Pa"),
                umag  = qty(m -> [sqrt(sum(abs2, u)) for u ∈ m.s.u], "displacement magnitude", "m"; label=L"u", cblim=(0.0, 0.25)),
                vmag  = qty(m -> [sqrt(sum(abs2, v)) for v ∈ m.s.v], "velocity magnitude", "m/s"; label=L"v", cblim=(0.0, 5.0)),
                c0    = qty(m -> [c.c₀ for c ∈ m.s.cmp], "initial cohesion", "Pa";
                            label=L"c_{0}", display=("kPa", 1e-3), cb=:vik, cblim=(10.0, 30.0)),
                phi0  = qty(m -> [hasproperty(c, :ϕ₀) ? c.ϕ₀ : NaN for c ∈ m.s.cmp], "initial friction angle", "rad";   # NaN: von Mises
                            label=L"\phi_{0}", display=("deg", 180/π), cblim=(20.0, 40.0)),
            ),
            fluid = (;),                                         # until the fluid phase has fields
            thermal = (; T = qty(m -> m.t.T, "temperature", "K"; label=L"T", cb=:thermal, cblim=(20.0, 293.15))),
        ),
        mesh = (;
            solid = (;
                v = qty(g -> columns(g.s.v), "velocity", "m/s"),
                a = qty(g -> columns(g.s.a), "acceleration", "m/s^2"),
                m = qty(g -> g.s.m, "mass", "kg"),
            ),
            thermal = (; T = qty(g -> g.t.T, "temperature", "K")),
        ),
    )
end

# the entry at a path ("point/solid/p"), or a group of them ("point/solid"); nothing if there is none
function quantity(path::AbstractString)
    t = quantities()
    for k ∈ Symbol.(split(path, "/"))
        (t isa NamedTuple && haskey(t, k)) || return nothing
        t = getfield(t, k)
    end
    return t
end
# the plottable point quantities, key => entry, in order: what `plot.what` picks from
function plottable()
    return [string(k) => e for phase ∈ quantities().point for (k, e) ∈ pairs(phase) if !isnothing(e.label)]
end
# a plottable quantity by key, or nothing
function plot_entry(key)
    for (k, e) ∈ plottable()
        k == string(key) && return e
    end
    return nothing
end
# a value of `plot.what` or of an `output.what` list (option path ["output", "what", "point", "solid"]) as menus
# show it: its quantity's readable name, otherwise the value itself
function shown_name(path, v)
    e = path[1] == "plot" ? plot_entry(v) : path[1] == "output" ? quantity(join([path[3:end]; string(v)], "/")) : nothing
    return isnothing(e) ? string(v) : e.name
end
# the `output.what` options of `get_option()`: per object and phase, the keys of its quantities
function output_what_options()
    return NamedTuple(object => NamedTuple(phase => ("Select $phase $object fields", String[string(k) for k ∈ keys(qs)])
                                           for (phase, qs) ∈ pairs(phases))
                      for (object, phases) ∈ pairs(quantities()))
end
