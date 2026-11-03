export draw_case, zone_at, zone_index, property, points, partition, write_geo, save_layout   # EasyDraw's, on a case too

"""
    draw_case(; L=[10, 5], nel=[20, 10], ni=2) -> Case or nothing
    draw_case(case::Case)
    draw_case(file)

Set up a case in EasyDraw's window (`EasyDraw.draw`) on the case's layout and settings, in four steps, the geometry
first, then the physics: **Layout** (optionally start from a case file; the domain: draw on a domain, dragging the
handle off its corner or typing `Lx`, `Lz`, or on a picture, scaled by a scale bar measured on it, which the domain then
follows; the element size; a width along y, which makes the case 3D; points per element; and the zones, drawn as
polygons, rectangles, triangles or circles), **Physics** (each zone's optional cohesion, residual cohesion, friction
angle, temperature, Young's modulus, Poisson's ratio and density; the defaults, the case's material; a boundary
condition per side), **Options** (another catalog file can be loaded there; solver options from drop-down lists, any
value for numeric ones, ticked to override the defaults; the plotted variables) and **Create** (the run settings:
workflow, times or load steps; `check_case` and the estimated time step).
What it offers is `get_catalog()`. **Create…**, at the end of the Create step, asks for the case's folder (chosen or made
in the system's folder dialog): the case is written there as `<name>_case.toml`, its runs go there (`run.folder`:
`execute!(setup_simulation(case))`), and the window returns the `Case` (errors block it, warnings ask for a second
click). **Cancel**, closing the window or Ctrl-C returns `nothing`. It opens on the Layout step. It draws zones in x–z
(polygons, blocks, balls as circles).
"""
function draw_case(case::Case; kwargs...)
    r = EasyDraw.draw(case.layout, case.settings; catalog=get_catalog(), case_hooks()..., kwargs...)
    return isnothing(r) ? nothing : Case(r.layout, r.settings)
end
function draw_case(file::AbstractString; kwargs...)
    return draw_case(load_case(file); kwargs...)
end
function draw_case(; L=[10.0, 5.0], nel=[20, 10], ni::Integer=2, kwargs...)
    return draw_case(Case(Layout(; L, nel, ni), case_settings(; run=(; folder="drawn"))); kwargs...)
end

"""
    case_hooks() -> NamedTuple

What the drawing window calls back for a case: its checks (`check_case`), what **Create…** writes (the case as
`<name>_case.toml` in the chosen folder, which becomes its run folder), what "Start from a file…" reads (a case file),
and its title.
"""
function case_hooks()
    check = (l, s) -> (c = Case(l, s); check_case(c, Solver(c)))
    # Create…: the case's runs go to the chosen folder. EasyDraw returns these same settings, so the case drawn has it too.
    save  = (l, s, dir) -> (s.run["folder"] = dir; save_case(joinpath(dir, "$(basename(dir))_case.toml"), Case(l, s)))
    load  = file -> (c = load_case(file); (c.layout, c.settings))
    return (; title="ElastoPlasm — draw a case", check, save, load)
end
