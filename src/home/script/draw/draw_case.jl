export draw_case, assign_case, zone_at, zone_index, property, points, partition, write_geo, save_layout   # EasyDraw's, on a case too

"""
    draw_case(; L=[10, 5], nel=[20, 10], ni=2) -> Case or nothing
    draw_case(case::Case)
    draw_case(file)

Set up a case in EasyDraw's window (`EasyDraw.draw`), on its four pages, the geometry first, then the physics:
**Layout** (the domain: drawn on a domain, dragging the handle off its corner or typing `Lx`, `Lz`, or on a picture,
scaled by a scale bar measured on it, which the domain then follows; the element size; a width along y, which makes
the case 3D; points per element; and the zones, drawn as polygons, rectangles, triangles or circles), **Physics**
(each zone's optional cohesion, residual cohesion, friction angle, temperature, Young's modulus, Poisson's ratio and
density; the defaults, the case's material; a boundary condition per side), **Options** (another catalog file can be
loaded there; solver options from drop-down lists, any value for numeric ones, ticked to override the defaults; the
plotted variables) and **Create** (the run settings: workflow, times or load steps; `check_case` and the estimated
time step). "Start from a file…" opens a case file.
What it offers is `get_catalog()`. **Create…**, at the end of the Create page, asks for the case's folder (chosen or
made in the system's folder dialog): the case is written there as `<name>_case.toml`, its runs go there (`run.folder`:
`execute!(setup_simulation(case))`), and the window returns the `Case` (errors block it, warnings ask for a second
click). **Cancel**, closing the window or Ctrl-C returns `nothing`. It draws zones in x–z (polygons, blocks, balls as
circles). To give the same geometry other values, see `assign_case`.
"""
function draw_case(case::Case; kwargs...)
    return EasyDraw.draw(case; catalog=get_catalog(), case_hooks()..., kwargs...)
end
function draw_case(file::AbstractString; kwargs...)
    return draw_case(load_case(file); kwargs...)
end
function draw_case(; L=[10.0, 5.0], nel=[20, 10], ni::Integer=2, kwargs...)
    return draw_case(Case(Layout(; L, nel, ni), case_settings(; run=(; folder="drawn"))); kwargs...)
end

"""
    assign_case(case::Case) -> Case or nothing
    assign_case(file)

Give a case's geometry other values in EasyDraw's window (`EasyDraw.assign`): its **Physics**, **Options** and
**Create** pages, the layout shown but not edited, so that one geometry gets several settings. Create… writes the case
as `draw_case` does; "Start from a file…" takes a case file's settings (its zones must be the layout's).
"""
function assign_case(case::Case; kwargs...)
    h = case_hooks()
    load = file -> Case(case.layout, load_case(file).settings)
    return EasyDraw.assign(case; catalog=get_catalog(), h..., load, kwargs...)
end
function assign_case(file::AbstractString; kwargs...)
    return assign_case(load_case(file); kwargs...)
end

"""
    case_hooks() -> NamedTuple

What the drawing window calls back for a case: its checks (`check_case` with its solver), what **Create…** writes (the
case, its run folder set to the chosen folder, as `<name>_case.toml` there; the window returns that case), what "Start
from a file…" reads (a case file), and its title.
"""
function case_hooks()
    check = c -> check_case(c, Solver(c))
    save  = function (c, dir)
        run = merge(c.settings.run, Dict{String,Any}("folder" => dir))             # its runs go to the chosen folder
        s = c.settings
        written = Case(c.layout, Settings(; zones=s.zones, defaults=s.defaults, sides=s.sides, options=s.options, run))
        save_case(joinpath(dir, "$(basename(dir))_case.toml"), written)
        return written
    end
    return (; title="ElastoPlasm — draw a case", check, save, load=load_case)
end
