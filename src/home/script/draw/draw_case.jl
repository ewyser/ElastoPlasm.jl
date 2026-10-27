export draw_case

"""
    draw_case(; L=[10, 5], nel=[20, 10], ni=2) -> Case or nothing
    draw_case(case::Case)
    draw_case(file)

Set up a case in a window, in four steps: **Domain** (optionally start from a case file; draw on a domain, dragging
its corner or typing `Lx`, `Lz`, or on a picture, scaled by a scale bar measured on it, which the domain then
follows; the element size; a width to extrude along y into 3D; points per element; a boundary condition per side),
**Zones** (draw the material as polygons, each with optional cohesion, residual cohesion, friction
angle, temperature, Young's modulus, Poisson's ratio and density), **Options** (solver options from
drop-down lists, any value for numeric ones, ticked to override the defaults; the plotted variables)
and **Run** (workflow, times or load steps; `check_case` and the estimated time step). **Create…**, at the end of
the Run step, asks for the case's folder (chosen or made in the system's folder dialog): the case is written there
as TOML, its runs go there (`run.folder`: `execute!(setup_simulation(case))`), and
the window returns the `Case` (errors block it, warnings ask for a second click). **Cancel**, closing the window or
Ctrl-C returns `nothing`. Without arguments the
window opens on the domain step; with a case (or its file), on the zones step. It edits 2D (x–z)
cases with polygon, block or ball (circle) zones.

The window comes from the QML extension (`DrawQMLExt`, Qt Quick): `using QML`. Install QML in your
default environment with Qt pinned to the version Plots/GR already load (a QML built for another Qt fails
to load): without the window, `draw_case()`'s error gives the command with that version.
"""
draw_case(case::Case) = drawing_window(case; step=2)
draw_case(file::AbstractString) = draw_case(load_case(file))
draw_case(; L=[10.0, 5.0], nel=[20, 10], ni::Integer=2) = drawing_window(Case(; L, nel, ni, run=(; folder="drawn")); step=1)

"""
    drawing_window(case; step) -> Case or nothing

The window `draw_case` opens on `case` at `step` (1 domain … 4 run), defined by the QML extension; without
it, an error saying how to get it (`drawing_window_help`).
"""
drawing_window(case; step) = error(drawing_window_help())

# the Qt this session already loads (Plots/GR load Qt6Base_jll), which QML must be built for: "6.8.2", or nothing
function session_qt()
    m = [v for (k, v) ∈ Base.loaded_modules if k.name == "Qt6Base_jll"]
    isempty(m) && return nothing
    v = pkgversion(only(m))
    return "$(v.major).$(v.minor).$(v.patch)"
end

"""
    drawing_window_help(; installed, qt) -> String

What to do when the drawing window is missing: install QML, or load it, with every Qt library pinned to
`qt`, the version this session loads already. A QML built for another Qt fails at `using QML` with a missing
`_qt_version_tag_…` symbol, before ElastoPlasm can say anything, hence the explanation here.
"""
function drawing_window_help(; installed::Bool=!isnothing(Base.find_package("QML")), qt=session_qt())
    v    = something(qt, "<the Qt version of Qt6Base_jll in the project's Manifest>")
    pins = join(("$(p)@$(v)" for p ∈ ("Qt6Base_jll", "Qt6Declarative_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll", "Qt6Wayland_jll")), " ")
    fix  = "] activate; add QML $pins"
    why  = "a missing `_qt_version_tag` symbol, QML was built for another Qt than the $v this session loads through Plots/GR"
    msg  = installed ? "QML is installed but not loaded: run `using QML`. If it fails with $why; pin them together: $fix." :
                       "install QML in your default environment, with Qt pinned to the $v this session loads through Plots/GR " *
                       "(otherwise `using QML` fails with a missing `_qt_version_tag` symbol): $fix, then `using QML`."
    return "draw_case needs the QML extension: " * msg * " Without the window, build a Case in code or write it as TOML."
end
