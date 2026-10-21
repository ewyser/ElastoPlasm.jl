# What the QML view reads: the canvas and the side panel as JSON strings, in the canvas's pixels.
# Independent of Qt (test_draw.jl checks it without QML).

# ---- JSON for QML (the canvas and the panel read their contents as JSON strings)

# QML counts from 0, Julia from 1: every index crosses the boundary here, out with `qml_index` (no item → -1)
# and back with `julia_index`; the QML file does no index arithmetic
qml_index(i) = isnothing(i) ? -1 : i - 1
fileurl(f) = "file://" * replace(abspath(f), "%" => "%25", " " => "%20", "#" => "%23")   # a local file, for QML
julia_index(i) = Int(i) + 1

json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\t" => "\\t") * "\""
json(x::Symbol)  = json(string(x))
json(x::Bool)    = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::Real)    = isfinite(x) ? string(round(Float64(x); sigdigits=10)) : "null"
json(::Nothing)  = "null"
json(x::Union{AbstractVector,Tuple}) = "[" * join(json.(x), ",") * "]"
json(x::NamedTuple)   = "{" * join(("$(json(string(k))):$(json(v))" for (k, v) ∈ pairs(x)), ",") * "}"
json(x::AbstractDict) = "{" * join(("$(json(string(k))):$(json(v))" for (k, v) ∈ x), ",") * "}"
struct RawJSON; s::String; end                              # already JSON (the cached points)
json(x::RawJSON) = x.s

# ---- pixels ↔ data: the visible area `view`, equal aspect, margins for the axis labels

const PAD = (; left=48.0, right=12.0, top=28.0, bottom=36.0)
function transform(e::Editor, w, h)
    x0, x1, z0, z1 = e.view
    aw, ah = w - PAD.left - PAD.right, h - PAD.top - PAD.bottom
    s  = min(aw/(x1 - x0), ah/(z1 - z0))
    ox = PAD.left + (aw - s*(x1 - x0))/2
    oy = PAD.bottom + (ah - s*(z1 - z0))/2
    topx(x, z) = (ox + (x - x0)*s, h - (oy + (z - z0)*s))
    todata(px, py) = (x0 + (px - ox)/s, z0 + (h - py - oy)/s)
    return (; topx, todata, s, box=(ox, h - oy - s*(z1 - z0), s*(x1 - x0), s*(z1 - z0)))
end
nice_tick(r) = (b = 10.0^floor(log10(r/5)); m = r/5/b; b*(m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10))

# the snap cursor: green on the first vertex of the zone being drawn (a click there closes it), red while a zone
# is drawn, grey otherwise; its shape says what it snapped to (a square: a vertex, a diamond: an edge, a ring:
# the lattice, or free)
function cursor_json(e::Editor, P)
    isnothing(e.cursor) && return nothing
    q, kind = e.cursor
    color = kind == :close ? "#2a9d6f" : isnothing(e.drawing) ? "#9ca3af" : "#dc2626"
    shape = kind == :vertex ? "square" : kind == :edge ? "diamond" : "ring"
    return (; xy=P(q...), r=kind == :close ? 9 : 6, color, shape)
end
"""
    background_json(editor, w, h) -> String
    overlay_json(editor, w, h) -> String

What the canvas draws, in pixels of a `w`×`h` canvas, in two layers. The background: the axes box,
ticks and title, the grid lines, the domain and its sides (fixed red, roller blue, free grey), the
zones' fills and the points to be created (over the fills, so each shows its zone's colour). The
overlay, on top: the zones' outlines and vertices, the polygon being drawn, the snap cursor and its
coordinates, and the corner handle (domain step). A mouse move repaints only the overlay, so the
points (tens of thousands on a fine mesh) are not parsed and drawn again; while vertices are dragged,
the fills and points follow on release (when the background repaints).
"""
# each zone's outline in pixels and its handles (a polygon's vertices; a circle's centre and radius handle)
circle_px(P, c, r; n=96, θ=0.0) = [P(c[1] + r*cos(θ + 2π*j/n), c[2] + r*sin(θ + 2π*j/n)) for j ∈ 0:n-1]
zone_paths(e::Editor, P) = [iscircle(rg) ?
    (; pts=circle_px(P, rg.centre, rg.radius), handles=[P(rg.centre...), P(radius_handle(rg)...)], color=COLORS[mod1(k, length(COLORS))], selected=k == e.selected) :
    (; pts=[P(v...) for v ∈ rg.vertices], handles=[P(v...) for v ∈ rg.vertices], color=COLORS[mod1(k, length(COLORS))], selected=k == e.selected)
    for (k, rg) ∈ enumerate(e.regions)]
# the zone being drawn: the polygon so far and the cursor (a triangle closed back to its first corner); a rectangle
# from its first corner to the cursor; a circle at the cursor's distance from its centre
function sketch_px(e::Editor, P)
    d = e.drawing
    (isnothing(d) || isempty(d)) && return []
    c = isnothing(e.cursor) ? nothing : e.cursor[1]
    if e.shape == :rectangle
        a = d[1]
        return isnothing(c) ? [P(a...)] : [P(a...), P(c[1], a[2]), P(c...), P(a[1], c[2]), P(a...)]
    end
    if e.shape != :circle
        line = [[P(v...) for v ∈ d]; isnothing(c) ? [] : [P(c...)]]
        return e.shape == :triangle && length(d) == 2 ? [line; [P(d[1]...)]] : line
    end
    isnothing(c) && return [P(d[1]...)]
    ring = circle_px(P, d[1], hypot(c[1] - d[1][1], c[2] - d[1][2]); θ=atan(c[2] - d[1][2], c[1] - d[1][1]))
    return [ring; [ring[1]]; [P(d[1]...)]]                    # from the cursor round, then the radius to the centre
end
function background_json(e::Editor, w, h)
    t = transform(e, w, h)
    P(x, z) = collect(t.topx(x, z))
    L, nel, ni = e.L, e.nel, e.ni
    lines(n, h_, a) = [a == 1 ? [P(i*h_, 0); P(i*h_, L[2])] : [P(0, i*h_); P(L[1], i*h_)] for i ∈ 0:n]
    sub = [lines(nel[1]*ni, L[1]/(nel[1]*ni), 1); lines(nel[2]*ni, L[2]/(nel[2]*ni), 2)]
    el  = [lines(nel[1], L[1]/nel[1], 1); lines(nel[2], L[2]/nel[2], 2)]
    ends = Dict("left" => (0, 0, 0, L[2]), "right" => (L[1], 0, L[1], L[2]), "bottom" => (0, 0, L[1], 0), "top" => (0, L[2], L[1], L[2]))
    sides = [(; seg=[P(ends[n][1:2]...); P(ends[n][3:4]...)],
                color=e.bcs[a, k] == :fixed ? "#b22222" : e.bcs[a, k] == :roller ? "#4682b4" : "#cccccc") for (n, a, k) ∈ SIDES]
    fills = zone_paths(e, P)
    pts = e.points ? points_json(e, t, w, h) : RawJSON("[]")
    b = e.backdrop                                           # the picture: its file URL and pixel rectangle
    backdrop = isnothing(b) || b.npx[1] == 0 ? nothing :
               (; url=fileurl(b.file), rect=[P(b.x0, b.z0 + height(b)); b.width*t.s; height(b)*t.s], opacity=b.opacity)
    x0, x1, z0, z1 = e.view
    tx, tz = nice_tick(x1 - x0), nice_tick(z1 - z0)
    ticks = (; x=[(; px=P(v, z0)[1], label=string(round(v; sigdigits=6))) for v ∈ ceil(x0/tx)*tx:tx:x1],
               z=[(; py=P(x0, v)[2], label=string(round(v; sigdigits=6))) for v ∈ ceil(z0/tz)*tz:tz:z1])
    title = "$(L[1]) × $(L[2]) m, $(nel[1]) × $(nel[2]) elements, $(ni)×$(ni) points per element"
    return json((; box=collect(t.box), sub, el, domain=[P(0, 0); P(L[1], L[2])], sides, fills, points=pts, backdrop, ticks, title))
end
function overlay_json(e::Editor, w, h)
    t = transform(e, w, h)
    P(x, z) = collect(t.topx(x, z))
    zones = zone_paths(e, P)
    sketch = sketch_px(e, P)
    picked = isnothing(e.picked) || isempty(shared(e.regions, e.picked)) ? nothing : P(e.picked...)
    return json((; box=collect(t.box), zones, sketch, cursor=cursor_json(e, P), handle=e.step == 1 ? P(e.L...) : nothing, picked,
                   readout=e.readout, calib=isnothing(e.calib) ? [] : [(; xy=P(p...), n=j) for (j, p) ∈ enumerate(e.calib)]))
end

# the points to be created, per zone, as x1, y1, x2, y2, … in pixels (0.1 px): cached until the points, the
# visible area or the canvas size change, since a fine mesh has tens of thousands of them and the cursor
# repaints the canvas at every move
function points_json(e::Editor, t, w, h)
    p = preview_points(e)
    return cached!(e, :points_json, hash((objectid(p), e.view, w, h))) do
        x, z, r = p
        zones = map(unique(r)) do k
            xy = join((string(round(v; digits=1)) for i ∈ findall(==(k), r) for v ∈ t.topx(x[i], z[i])), ",")
            "{\"color\":$(json(COLORS[mod1(k, length(COLORS))])),\"xy\":[$xy]}"
        end
        RawJSON("[" * join(zones, ",") * "]")
    end
end

"""
    panel_json(editor) -> String

What the side panel shows: the step, the status line, and the fields of each step (domain, zones,
the rows of the current options page, run settings and checks).
"""
function panel_json(e::Editor)
    defaults = Dict(string(k) => getfield(e.base.material, k) for k ∈ keys(MATERIAL))
    vm = getoption(e.options, "material", ["plastic"]) == "VM"
    i  = e.selected
    domain = (; Lx=string(e.L[1]), Lz=string(e.L[2]), h=string(round(e.h; sigdigits=8)), width=isnothing(e.width) ? "" : string(e.width),
                nel="→ $(e.nel[1]) × $(e.nel[2]) elements of $(round(e.L[1]/e.nel[1]; sigdigits=4)) × $(round(e.L[2]/e.nel[2]; sigdigits=4)) m" *
                    (isnothing(e.width) ? "" : ", × $(e.nely) along y"),
                ni=qml_index(e.ni), nis=["$(n)×$(n)" for n ∈ 1:max(8, e.ni)], kinds=string.(collect(BOUNDARY_KINDS)),
                bcs=[(; side=n, index=qml_index(findfirst(==(e.bcs[a, k]), collect(BOUNDARY_KINDS)))) for (n, a, k) ∈ SIDES])
    x, _, r = preview_points(e)
    zones = (; names=["$k  $(rg.name)" for (k, rg) ∈ enumerate(e.regions)], selected=qml_index(i), name=i == 0 ? "" : e.regions[i].name,
               props=[(; key, label, default="default $(round(defaults[key]*scale; sigdigits=6))",
                         value=(i == 0 || !haskey(e.regions[i].props, key)) ? "" : string(round(e.regions[i].props[key]*scale; sigdigits=8)))
                      for (key, label, scale) ∈ PROPS if !(vm && key == "phi")],
               circle=i == 0 || !iscircle(e.regions[i]) ? nothing :
                      (; x=string(round(e.regions[i].centre[1]; sigdigits=8)), z=string(round(e.regions[i].centre[2]; sigdigits=8)),
                         r=string(round(e.regions[i].radius; sigdigits=8))),
               vertex=isnothing(e.picked) || isempty(shared(e.regions, e.picked)) ? nothing :
                      (; x=string(round(e.picked[1]; sigdigits=8)), z=string(round(e.picked[2]; sigdigits=8))),
               grid=e.grid, shapes=e.shapes, points=e.points, drawing=!isnothing(e.drawing), shape=string(e.shape),
               summary="$(length(x)) material points" * (isempty(e.regions) ? "" : ": ") *
                       join(("$(rg.name) $(count(==(k), r))" for (k, rg) ∈ enumerate(e.regions)), ", "))
    options = (; sections=option_sections(), section=qml_index(findfirst(==(e.section), option_sections())), rows=option_rows_json(e))
    c  = checks(e)
    qs = ElastoPlasm.quasistatic(e.run["workflow"])
    fields = qs ? [("lstps", "load steps")] : [("te", "elastic stage te [s]"), ("tg", "gravity ramp tg [s]"), ("tep", "plastic stage tep [s]")]
    run = (; workflows=collect(WORKFLOWS), workflow=qml_index(findfirst(==(e.run["workflow"]), collect(WORKFLOWS))),
             fields=[(; key, label, value=string(e.run[key])) for (key, label) ∈ [fields; ("fid", "folder in dump/")]],
             errors=join("✗ " .* c.errors, "\n"), warnings=join("! " .* c.warnings, "\n"), info=join(c.info, "\n"))
    b = e.backdrop
    backdrop = isnothing(b) ? nothing : (; name=basename(b.file), url=fileurl(b.file), loaded=b.npx[1] > 0,
                                           x0=fmt(round(b.x0; sigdigits=8)), z0=fmt(round(b.z0; sigdigits=8)),
                                           width=fmt(round(b.width; sigdigits=8)), opacity=fmt(b.opacity))
    calib = isnothing(e.calib) ? nothing :                   # a scale bar being measured: its ends clicked, its length
            (; clicked=length(e.calib), length=isnothing(e.calib_len) ? "" : fmt(e.calib_len),
               shown=length(e.calib) == 2 ? fmt(round(hypot((e.calib[2] .- e.calib[1])...); sigdigits=6)) : "")
    return json((; tab=qml_index(e.step), status=e.status, finished=e.finished, domain, zones, options, run, backdrop, calib))
end
# the rows of the current options page: a drop-down, a numeric field (with presets) or a list of values
function option_rows_json(e::Editor)
    rows = []
    for (sec, path, vals, def) ∈ filter(r -> option_page(r[1], r[2]) == e.section, option_rows())
        cur = getoption(e.options, sec, path)
        base = (; key=option_key(sec, path), label=join(path, "."), checked=!isnothing(cur), values=option_label.(vals))
        push!(rows, free_numeric(def) ? (; base..., kind="numeric", text=fmt(something(cur, def)), index=-1) :
                                        (; base..., kind="menu", text="", index=qml_index(something(findfirst(v -> option_label(v) == option_label(something(cur, def)), vals), 1))))
    end
    for (sec, path, vals, def) ∈ filter(r -> option_page(r[1], r[2]) == e.section, option_lists())
        cur  = getoption(e.options, sec, path)
        vals = list_values(sec, vals)
        on   = option_label.(something(cur, def))
        push!(rows, (; key=option_key(sec, path), label=join(path, "."), checked=!isnothing(cur), values=option_label.(vals),
                       kind="list", text="", index=-1, picks=[option_label(v) ∈ on for v ∈ vals]))
    end
    return rows
end
