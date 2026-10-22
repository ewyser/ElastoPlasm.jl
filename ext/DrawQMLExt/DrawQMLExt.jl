# QML (Qt Quick) window for draw_case: a thin view of an `Editor` (editor.jl). Install QML in your
# default environment (] activate; add QML) and `using QML` next to ElastoPlasm to enable it.
module DrawQMLExt

using ElastoPlasm, QML
include("editing.jl")                                        # the window without Qt: Editor and the JSON it shows

const EDITOR = Ref{Editor}()                                 # the editor the open window shows
const QMLFILE = joinpath(@__DIR__, "DrawCase.qml")

# ---- what QML calls (`Julia.<name>(…)`), on the open window's editor: QML's strings arrive as Julia strings, its
# 0-based indices go through `julia_index`; an action returns nothing (QML needs no editor back)

plain(x) = x isa AbstractString ? String(x) : x
act(f) = (args...) -> (f(EDITOR[], plain.(args)...); nothing)
todata(e, px, py, w, h) = transform(e, w, h).todata(px, py)
const CALLS = (;
    background = (w, h) -> background_json(EDITOR[], Float64(w), Float64(h)),
    overlay    = (w, h) -> overlay_json(EDITOR[], Float64(w), Float64(h)),
    panel      = () -> panel_json(EDITOR[]),
    press      = act((e, px, py, right, w, h) -> press!(e, todata(e, px, py, w, h)...; right=Bool(right))),
    move       = act((e, px, py, w, h) -> move!(e, todata(e, px, py, w, h)...)),
    release    = act(release!),
    leave      = act(leave!),
    wheel      = act((e, px, py, delta, w, h) -> zoom!(e, todata(e, px, py, w, h)..., 0.9^(Float64(delta)/120))),   # about the pointer
    key        = act((e, name, cmd, shift) -> key!(e, name; cmd=Bool(cmd), shift=Bool(shift))),
    set_step   = act((e, i) -> set_step!(e, julia_index(i))),                     # the tab's index
    set_field  = act(set_domain_field!),
    set_ni     = act((e, i) -> e.ni = julia_index(i)),
    set_bc     = act((e, side, j) -> set_bc!(e, side, BOUNDARY_KINDS[julia_index(j)])),
    new_zone   = act((e, shape) -> new_zone!(e, Symbol(shape))),                  # "polygon", "rectangle", "triangle", "circle"
    set_circle = act(set_circle!),
    select_zone = act((e, i) -> select_zone!(e, julia_index(i))),
    rename_zone = act(rename_zone!),
    set_prop   = act(set_prop!),
    swap_zone  = act((e, d) -> swap_zone!(e, Int(d))),
    delete_zone = act(delete_zone!),
    set_picked = act(set_picked!),
    remove_picked = act(remove_picked!),
    undo       = act(undo!),
    redo       = act(redo!),
    toggle     = act((e, name, on) -> setfield!(e, Symbol(name), Bool(on))),     # the snapping and points check boxes
    set_section = act((e, i) -> e.section = SECTIONS[julia_index(i)]),
    set_option = act((e, k, j) -> set_option!(e, k, julia_index(j))),
    type_option = act(type_option!),
    tick_option = act((e, k, on) -> tick_option!(e, k, Bool(on))),
    tick_value = act((e, k, i, on) -> tick_value!(e, k, julia_index(i), Bool(on))),
    set_workflow = act((e, j) -> set_run!(e, "workflow", WORKFLOWS[julia_index(j)])),
    set_run    = act(set_run!),
    ready      = () -> ready!(EDITOR[]),
    create     = url -> create!(EDITOR[], localpath(String(url))),
    cancel     = act(cancel!),
    load       = act((e, url) -> load!(e, localpath(url))),
    set_backdrop = act((e, url) -> set_backdrop!(e, localpath(url))),
    remove_backdrop = act(remove_backdrop!),
    backdrop_size = act(backdrop_size!),
    set_backdrop_field = act(set_backdrop_field!),
    scale_bar  = act(scale_bar!),
    set_scale_length = act(set_scale_length!),
    apply_scale = act(apply_scale!),
    fit_domain = act(fit_domain!),
)

# registered for each window: QML.jl forgets them when a window's application ends, so registering once
# (in __init__) left a second window of the same session with none
register_functions() = foreach(((name, f),) -> QML.qmlfunction(string(name), f), pairs(CALLS))

ElastoPlasm.drawing_window(case::Case; step::Integer=2) = qml_window(case; step)

"""
    qml_window(case; step=2, shots=nothing, startdir=pwd()) -> Case or nothing

Open the QML window on `case` (blocking until Create…, Cancel or closing it) and return the edited case,
or `nothing`. With `shots = (dir, steps)` the window instead saves a PNG of each listed step into `dir`
and closes (for tests and documentation; run it with `QT_QPA_PLATFORM=offscreen`). Its file and folder
dialogs start in `startdir`: the current directory, not `dump/`, the solver's scratch folder
(`clean_dump()` empties it).
"""
function qml_window(case::Case; step::Integer=2, shots=nothing, startdir=pwd())
    EDITOR[] = Editor(case; step)
    isnothing(shots) || (ENV["QT_QUICK_CONTROLS_STYLE"] = get(ENV, "QT_QUICK_CONTROLS_STYLE", "Fusion"))
    register_functions()
    steps = isnothing(shots) ? Int[] : collect(Int, shots[2])
    loadqml(QMLFILE; shotdir=isnothing(shots) ? "" : shots[1], shotsteps=steps, shottabs=qml_index.(steps),   # names, tabs
            startdir=fileurl(startdir))
    @info "Drawing window open: domain, zones, options and run, then click Create… (Cancel or closing the window aborts)"
    exec()
    return EDITOR[].result
end

end
