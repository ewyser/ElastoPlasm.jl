# QML (Qt Quick) window for draw_case: a thin view of an `Editor` (editor.jl). Install QML in your
# default environment (] activate; add QML) and `using QML` next to ElastoPlasm to enable it.
module DrawQMLExt

using ElastoPlasm, QML
include("editing.jl")                                        # the window without Qt: Editor and the JSON it shows

const EDITOR = Ref{Editor}()                                 # the editor the open window shows
const QMLFILE = joinpath(@__DIR__, "DrawCase.qml")

# ---- what QML calls (on the open window's editor)

ed() = EDITOR[]
todata(px, py, w, h) = transform(ed(), w, h).todata(px, py)
background(w, h) = background_json(ed(), Float64(w), Float64(h))
overlay(w, h)    = overlay_json(ed(), Float64(w), Float64(h))
panel()     = panel_json(ed())
press(px, py, right, w, h) = (press!(ed(), todata(px, py, w, h)...; right=Bool(right)); nothing)
move(px, py, w, h)         = (move!(ed(), todata(px, py, w, h)...); nothing)
release()                  = (release!(ed()); nothing)
leave()                    = (leave!(ed()); nothing)
wheel(px, py, delta, w, h) = (zoom!(ed(), todata(px, py, w, h)..., 0.9^(Float64(delta)/120)); nothing)   # about the pointer
key(name, cmd, shift)      = (key!(ed(), String(name); cmd=Bool(cmd), shift=Bool(shift)); nothing)
set_step(i)                = (set_step!(ed(), julia_index(i)); nothing)   # the tab's index
function set_field(name, text)                               # the domain fields
    e, t = ed(), String(text)
    if name == "width"
        set_width!(e, t); set_domain!(e, e.L)
    else
        x = tryparse(Float64, strip(t)); (isnothing(x) || x ≤ 0) && return nothing
        name == "h"  ? set_domain!(e, e.L; h=x) : name == "Lx" ? set_domain!(e, [x, e.L[2]]) : set_domain!(e, [e.L[1], x])
        fit!(e)
    end
    return nothing
end
set_ni(i)          = (ed().ni = julia_index(i); nothing)
set_bc(side, j)    = (set_bc!(ed(), String(side), BOUNDARY_KINDS[julia_index(j)]); nothing)
new_zone(shape)    = (new_zone!(ed(), Symbol(String(shape))); nothing)   # "polygon", "rectangle", "triangle", "circle"
set_circle(k, s)   = (set_circle!(ed(), String(k), String(s)); nothing)
select_zone(i)     = (select_zone!(ed(), julia_index(i)); nothing)
rename_zone(s)     = (rename_zone!(ed(), String(s)); nothing)
set_prop(k, s)     = (set_prop!(ed(), String(k), String(s)); nothing)
swap_zone(d)       = (swap_zone!(ed(), Int(d)); nothing)
delete_zone()      = (delete_zone!(ed()); nothing)
set_picked(a, s)   = (set_picked!(ed(), String(a), String(s)); nothing)
remove_picked()    = (remove_picked!(ed()); nothing)
undo()             = (undo!(ed()); nothing)
redo()             = (redo!(ed()); nothing)
toggle(name, on)   = (setfield!(ed(), Symbol(name), Bool(on)); nothing)
set_section(i)     = (ed().section = option_sections()[julia_index(i)]; nothing)
set_option(k, j)   = (set_option!(ed(), String(k), julia_index(j)); nothing)
type_option(k, s)  = (type_option!(ed(), String(k), String(s)); nothing)
tick_option(k, on) = (tick_option!(ed(), String(k), Bool(on)); nothing)
tick_value(k, i, on) = (tick_value!(ed(), String(k), julia_index(i), Bool(on)); nothing)
set_workflow(j)    = (set_run!(ed(), "workflow", WORKFLOWS[julia_index(j)]); nothing)
set_run(k, s)      = (set_run!(ed(), String(k), String(s)); nothing)
done()             = done!(ed())
cancel()           = (ed().finished || cancel!(ed()); nothing)
localpath(url)     = (u = String(url); startswith(u, "file://") ? replace(u[8:end], "%20" => " ", "%23" => "#", "%25" => "%") : u)
save(url)          = (save!(ed(), localpath(url)); nothing)
load(url)          = (load!(ed(), localpath(url)); nothing)
set_backdrop(url)  = (set_backdrop!(ed(), localpath(url)); nothing)
remove_backdrop()  = (remove_backdrop!(ed()); nothing)
backdrop_size(w, h) = (backdrop_size!(ed(), w, h); nothing)
set_backdrop_field(k, s) = (set_backdrop_field!(ed(), String(k), String(s)); nothing)
scale_bar()        = (scale_bar!(ed()); nothing)
set_scale_length(s) = (set_scale_length!(ed(), String(s)); nothing)
apply_scale()      = (apply_scale!(ed()); nothing)
fit_domain()       = (fit_domain!(ed()); nothing)

# registered for each window: QML.jl forgets them when a window's application ends, so registering once
# (in __init__) left a second window of the same session with none
register_functions() = @qmlfunction background overlay panel press move release leave wheel key set_step set_field set_ni set_bc new_zone set_circle select_zone rename_zone set_prop swap_zone delete_zone set_picked remove_picked undo redo toggle set_section set_option type_option tick_option tick_value set_workflow set_run done cancel save load set_backdrop remove_backdrop backdrop_size set_backdrop_field scale_bar set_scale_length apply_scale fit_domain

ElastoPlasm.drawing_window(case::Case; step::Integer=2) = qml_window(case; step)

"""
    qml_window(case; step=2, shots=nothing, startdir=pwd()) -> Case or nothing

Open the QML window on `case` (blocking until Done, Cancel or closing it) and return the edited case,
or `nothing`. With `shots = (dir, steps)` the window instead saves a PNG of each listed step into `dir`
and closes (for tests and documentation; run it with `QT_QPA_PLATFORM=offscreen`). Save… and Load… start in
`startdir`: the current directory, not `dump/`, which is emptied whenever the package recompiles.
"""
function qml_window(case::Case; step::Integer=2, shots=nothing, startdir=pwd())
    EDITOR[] = Editor(case; step)
    isnothing(shots) || (ENV["QT_QUICK_CONTROLS_STYLE"] = get(ENV, "QT_QUICK_CONTROLS_STYLE", "Fusion"))
    register_functions()
    steps = isnothing(shots) ? Int[] : collect(Int, shots[2])
    loadqml(QMLFILE; shotdir=isnothing(shots) ? "" : shots[1], shotsteps=steps, shottabs=qml_index.(steps),   # names, tabs
            startdir="file://" * abspath(startdir))
    @info "Drawing window open: domain, zones, options and run, then click Done (Cancel or closing the window aborts)"
    exec()
    return EDITOR[].result
end

end
