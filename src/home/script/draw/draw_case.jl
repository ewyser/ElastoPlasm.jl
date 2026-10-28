export draw_case, zone_at, zone_index, property, points, partition, write_geo, save_layout

"""
    draw_case(; L=[10, 5], nel=[20, 10], ni=2) -> Case or nothing
    draw_case(case::Case)
    draw_case(file)

Set up a case in EasyDraw's window (`EasyDraw.draw`), in four steps: **Domain** (optionally start from a case file;
draw on a domain, dragging its corner or typing `Lx`, `Lz`, or on a picture, scaled by a scale bar measured on it,
which the domain then follows; the element size; a width to extrude along y into 3D; points per element; a boundary
condition per side), **Zones** (draw the material as polygons, rectangles, triangles or circles, each with optional
cohesion, residual cohesion, friction angle, temperature, Young's modulus, Poisson's ratio and density),
**Options** (solver options from drop-down lists, any value for numeric ones, ticked to override the defaults; the
plotted variables) and **Run** (workflow, times or load steps; `check_case` and the estimated time step). What it
offers is `get_catalog()`. **Create…**, at the end of the Run step, asks for the case's folder (chosen or made in the
system's folder dialog): the case is written there as `<name>_case.toml`, its runs go there (`run.folder`:
`execute!(setup_simulation(case))`), and the window returns the `Case` (errors block it, warnings ask for a second
click). **Cancel**, closing the window or Ctrl-C returns `nothing`. Without arguments the window opens on the domain
step; with a case (or its file), on the zones step. It edits 2D (x–z) cases with polygon, block or ball (circle) zones.
"""
function draw_case(case::Case; step::Integer=2, kwargs...)
    l = EasyDraw.draw(layout_of(case); catalog=get_catalog(), step, case_hooks()..., kwargs...)
    return isnothing(l) ? nothing : case_of(l)
end
draw_case(file::AbstractString; kwargs...) = draw_case(load_case(file); kwargs...)
draw_case(; L=[10.0, 5.0], nel=[20, 10], ni::Integer=2, kwargs...) = draw_case(Case(; L, nel, ni, run=(; folder="drawn")); step=1, kwargs...)

"""
    case_hooks() -> NamedTuple

What the drawing window calls back for a case: its checks (`check_case`), what **Create…** writes (the case as
`<name>_case.toml` in the chosen folder), what "Start from a file…" reads (a case file), and its title.
"""
function case_hooks()
    return (; title="ElastoPlasm — draw a case", check=l -> check_case(case_of(l)),
              save=(l, dir) -> save_case(joinpath(dir, "$(basename(dir))_case.toml"), case_of(l)),
              load=file -> layout_of(load_case(file)))
end

# the window's sides and the axis and end (lower, upper) of `bcs` each one is
const SIDE_AXES = ("left" => (1, 1), "right" => (1, 2), "bottom" => (2, 1), "top" => (2, 2))

"""
    layout_of(case) -> EasyDraw.Layout
    case_of(layout) -> Case

A 2D case as the drawing window edits it, and back: the domain, zones and points per element as they are; the
boundary conditions as the sides' kinds; the material as the zone constants' defaults; the solver overrides as its
options; the run settings, with the thermal flag (carried through, not shown); an absolute run folder as its folder,
which Create… sets.
"""
function layout_of(case::Case)
    length(case.L) == 2 || throw(ArgumentError("EasyDraw's layouts are 2D (x–z): a 3D case converts only when extruded along y (width, nely)"))
    run = Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(case.run))
    run["thermal"] = case.thermal
    return EasyDraw.Layout(; L=case.L, nel=case.nel, ni=case.ni, width=case.width, nely=case.nely, zones=case.zones,
                             sides=Dict(s => string(case.bcs[a, k]) for (s, (a, k)) ∈ SIDE_AXES),
                             defaults=Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(case.material)),
                             options=config_dict(case.solver), run, folder=isabspath(case.run.folder) ? case.run.folder : "")
end
function case_of(l::EasyDraw.Layout)
    run = copy(l.run)
    thermal = pop!(run, "thermal", false)
    isempty(l.folder) || (run["folder"] = l.folder)
    bcs = fill(:roller, 2, 2)
    foreach(((s, (a, k)),) -> bcs[a, k] = Symbol(get(l.sides, s, "roller")), SIDE_AXES)
    material = Dict{String,Any}(k => v for (k, v) ∈ l.defaults if haskey(MATERIAL, Symbol(k)))
    return Case(; L=l.L, nel=l.nel, ni=l.ni, width=l.width, nely=l.nely, zones=l.zones, bcs, material, solver=l.options,
                  thermal, run)
end

# ---- a case wherever EasyDraw takes a layout (`layout_of`; 2D, or extruded along y): its zone at a point (also 3D), on
# a grid, its points (those of the simulation, to the bit), its partition, its Gmsh script, its layout file (for
# EasyDraw alone; `load_case` reads it back). Exported with `using ElastoPlasm`.

function EasyDraw.zone_at(c::Case, x)
    return zone_at(c.zones, x)
end
function EasyDraw.zone_index(c::Case, xs, zs)
    return zone_index(layout_of(c), xs, zs)
end
function EasyDraw.property(c::Case, key, xs, zs; kwargs...)
    return property(layout_of(c), key, xs, zs; kwargs...)
end
function EasyDraw.points(c::Case)
    return points(layout_of(c))
end
function EasyDraw.partition(c::Case; kwargs...)
    return partition(layout_of(c); kwargs...)
end
function EasyDraw.write_geo(file::AbstractString, c::Case; kwargs...)
    return write_geo(file, layout_of(c); kwargs...)
end
function EasyDraw.save_layout(file::AbstractString, c::Case)
    return save_layout(file, layout_of(c))
end
