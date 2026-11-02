export draw_case, zone_at, zone_index, property, points, partition, write_geo, save_layout

"""
    draw_case(; L=[10, 5], nel=[20, 10], ni=2) -> Case or nothing
    draw_case(case::Case)
    draw_case(file)

Set up a case in EasyDraw's window (`EasyDraw.draw`) on the case's layout, in four steps: **Domain** (optionally start
from a case file; draw on a domain, dragging its corner or typing `Lx`, `Lz`, or on a picture, scaled by a scale bar
measured on it, which the domain then follows; the element size; a width along y, which makes the case 3D; points per
element; a boundary condition per side), **Zones** (draw the material as polygons, rectangles, triangles or circles,
each with optional cohesion, residual cohesion, friction angle, temperature, Young's modulus, Poisson's ratio and
density), **Options** (solver options from drop-down lists, any value for numeric ones, ticked to override the
defaults; the plotted variables) and **Run** (workflow, times or load steps; `check_case` and the estimated time step).
What it offers is `get_catalog()`. **Create…**, at the end of the Run step, asks for the case's folder (chosen or made
in the system's folder dialog): the case is written there as `<name>_case.toml`, its runs go there (`run.folder`:
`execute!(setup_simulation(case))`), and the window returns the `Case` (errors block it, warnings ask for a second
click). **Cancel**, closing the window or Ctrl-C returns `nothing`. Without arguments the window opens on the domain
step; with a case (or its file), on the zones step. It draws zones in x–z (polygons, blocks, balls as circles).
"""
function draw_case(case::Case; step::Integer=2, kwargs...)
    l = EasyDraw.draw(case.layout; catalog=get_catalog(), step, case_hooks()..., kwargs...)
    return isnothing(l) ? nothing : Case(l)
end
function draw_case(file::AbstractString; kwargs...)
    return draw_case(load_case(file); kwargs...)
end
function draw_case(; L=[10.0, 5.0], nel=[20, 10], ni::Integer=2, kwargs...)
    return draw_case(Case(; L, nel, ni, run=(; folder="drawn")); step=1, kwargs...)
end

"""
    case_hooks() -> NamedTuple

What the drawing window calls back for a case: its checks (`check_case`), what **Create…** writes (the case as
`<name>_case.toml` in the chosen folder, which becomes its run folder), what "Start from a file…" reads (a case file),
and its title.
"""
function case_hooks()
    check = l -> check_case(Case(l))
    # Create…: the case's runs go to the chosen folder. EasyDraw returns this same layout, so the case drawn has it too.
    save  = (l, dir) -> (l.run["folder"] = dir; save_case(joinpath(dir, "$(basename(dir))_case.toml"), Case(l)))
    load  = file -> load_case(file).layout
    return (; title="ElastoPlasm — draw a case", check, save, load)
end

# ---- a case wherever EasyDraw takes a layout, through its own: its zone at a point, on a grid, its points (those of the
# simulation, to the bit, in 2D), its partition, its Gmsh script, its file. Exported with `using ElastoPlasm`.

function EasyDraw.zone_at(c::Case, x)
    return zone_at(c.layout, x)
end
function EasyDraw.zone_index(c::Case, xs, zs)
    return zone_index(c.layout, xs, zs)
end
# a property's default is the case's material constant
function EasyDraw.property(c::Case, key, xs, zs; default=get(c.material, Symbol(key), NaN))
    return property(c.layout, key, xs, zs; default)
end
function EasyDraw.points(c::Case)
    return points(c.layout)
end
function EasyDraw.partition(c::Case; kwargs...)
    return partition(c.layout; kwargs...)
end
function EasyDraw.write_geo(file::AbstractString, c::Case; kwargs...)
    return write_geo(file, c.layout; kwargs...)
end
function EasyDraw.save_layout(file::AbstractString, c::Case)
    return save_layout(file, c.layout)
end
