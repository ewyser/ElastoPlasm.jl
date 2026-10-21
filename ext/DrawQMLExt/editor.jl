# The drawing window as plain Julia state: what is drawn, and what each user action does. A GUI
# (`DrawQMLExt`) only shows it and forwards events, so everything here is testable without one.
# Coordinates are data coordinates (x, z in m) unless named `px`.

"""
    Editor(case; step=2)

The state of a drawing window editing `case`: the zones (as polygons), the domain, the boundary
conditions, the solver overrides and run settings, undo/redo snapshots, the polygon being drawn,
the vertices being dragged, the current step (1 Domain, 2 Zones, 3 Options, 4 Run), the options page,
the selected zone, the snapping toggles, the visible area (`view`), the cursor, its coordinates
(`readout`, drawn on the canvas) and the status line. `cache` keeps the preview points and the checks
until what they depend on changes (`preview_points`, `checks`).
`result` is the `Case` after **Done** (`done!`), or `nothing` when cancelled; `finished` says whether
the window should close.
"""
mutable struct Editor
    base     :: Case                                         # its material and thermal flag are kept
    regions  :: Vector{DrawnRegion}
    L        :: Vector{Float64}
    nel      :: Vector{Int}
    h        :: Float64                                      # requested element size [m]
    ni       :: Int
    width    :: Union{Nothing,Float64}
    nely     :: Union{Nothing,Int}
    bcs      :: Matrix{Symbol}
    options  :: Dict{String,Any}
    run      :: Dict{String,Any}
    history  :: Vector{Vector{DrawnRegion}}
    future   :: Vector{Vector{DrawnRegion}}
    drawing  :: Union{Nothing,Vector{Vector{Float64}}}
    onedge   :: Vector{Vector{Float64}}                      # vertices of the drawing snapped onto edges
    drag     :: Vector{Tuple{Int,Int}}                       # (zone, vertex) moving together
    resizing :: Bool                                         # the domain corner is being dragged
    step     :: Int
    section  :: String
    selected :: Int
    grid     :: Bool
    shapes   :: Bool
    points   :: Bool
    view     :: NTuple{4,Float64}                            # x0, x1, z0, z1
    cursor   :: Union{Nothing,Tuple{Vector{Float64},Symbol}}
    status   :: String
    confirm  :: Bool                                         # warnings shown: the next Done goes through
    result   :: Union{Nothing,Case}
    finished :: Bool
    readout  :: String                                       # the cursor's coordinates
    cache    :: Dict{Symbol,Tuple{UInt,Any}}                 # name => (key, value)
end
function Editor(case::Case; step::Integer=2)
    length(case.L) == 2 || throw(ArgumentError("the window draws a 2D (x–z) case; a width (domain step) extrudes it along y"))
    e = Editor(case, drawn_regions(case), copy(case.L), copy(case.nel), minimum(case.L ./ case.nel), case.ni, case.width,
               case.nely, copy(case.bcs), config_dict(case.solver), Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(case.run)),
               [], [], nothing, [], [], false, Int(step), first(option_sections()), isempty(case.zones) ? 0 : 1,
               true, true, true, (0.0, 1.0, 0.0, 1.0), nothing, "", false, nothing, false,
               "", Dict{Symbol,Tuple{UInt,Any}}())
    fit!(e)
    return e
end

lattice(e::Editor) = (e.L[1]/(e.nel[1]*e.ni), e.L[2]/(e.nel[2]*e.ni))
tol(e::Editor)     = 0.015*max(e.view[2] - e.view[1], e.view[4] - e.view[3])   # pick and snap radius [m]: follows the zoom
"Zoom by `f` (< 1: in) about the point `(x, z)`, which stays where it is on screen."
function zoom!(e::Editor, x, z, f)
    x0, x1, z0, z1 = e.view
    e.view = (x + (x0 - x)*f, x + (x1 - x)*f, z + (z0 - z)*f, z + (z1 - z)*f)
    return e
end
fit!(e::Editor, f=e.step == 1 ? 1.4 : 1.0) = (e.view = (-0.02f*e.L[1], (f + 0.02)*e.L[1], -0.04f*e.L[2], (f + 0.04)*e.L[2]); e)
checkpoint!(e::Editor) = (push!(e.history, deepcopy(e.regions)); length(e.history) > 200 && popfirst!(e.history); empty!(e.future); e)
outside(e::Editor) = any(v -> !(0 ≤ v[1] ≤ e.L[1] && 0 ≤ v[2] ≤ e.L[2]), (v for rg ∈ e.regions for v ∈ rg.vertices))
snap(e::Editor, p; skip=()) = snap_point(p, e.regions, lattice(e), e.L, tol(e); grid=e.grid, shapes=e.shapes, skip)
# the zone being drawn: is `p` on its first vertex (a click there closes it, from 3 vertices on)?
near_first(e::Editor, p) = !isnothing(e.drawing) && !isempty(e.drawing) &&
                           hypot(p[1] - e.drawing[1][1], p[2] - e.drawing[1][2]) < tol(e)
closing(e::Editor, p) = near_first(e, p) && length(e.drawing) ≥ 3
near_vertex(e::Editor, p) = (for rg ∈ reverse(e.regions), v ∈ rg.vertices
                                 hypot(p[1] - v[1], p[2] - v[2]) < tol(e) && return v
                             end; nothing)

# ---- the domain

function set_domain!(e::Editor, L; h=e.h)
    e.h, e.L, e.nel = h, Float64.(collect(L)), nel_for(L, h)
    isnothing(e.width) || (e.nely = nel_for([e.width], h)[1])
    e.status = outside(e) ? "Some zone vertices lie outside the domain." : ""
    return e
end
function set_width!(e::Editor, text)
    x = tryparse(Float64, strip(text))
    e.width = isnothing(x) || x ≤ 0 ? nothing : x
    e.nely  = isnothing(e.width) ? nothing : nel_for([e.width], e.h)[1]
    return e
end
set_bc!(e::Editor, side, kind) = (s = findfirst(t -> t[1] == side, SIDES); e.bcs[SIDES[s][2], SIDES[s][3]] = Symbol(kind); e)
function set_step!(e::Editor, k)
    e.step, e.drawing, e.confirm, e.cursor = Int(k), nothing, false, nothing
    fit!(e)
    return e
end

# ---- mouse and keys (data coordinates)

function press!(e::Editor, x, z; right::Bool=false)
    p = [Float64(x), Float64(z)]
    if e.step == 1                                           # the domain: grab the corner handle
        !right && hypot(p[1] - e.L[1], p[2] - e.L[2]) < 2tol(e) && (e.resizing = true)
        return e
    end
    e.step == 2 || return e
    if right                                                 # delete a vertex
        v = isnothing(e.drawing) ? near_vertex(e, p) : nothing
        isnothing(v) && return e
        snapshot = deepcopy(e.regions)
        delete_vertex!(e.regions, v) ? (push!(e.history, snapshot); empty!(e.future)) : (e.status = "A zone keeps at least 3 vertices.")
        return e
    end
    if !isnothing(e.drawing)                                 # drawing: a vertex, or close on the first one
        d = e.drawing
        if closing(e, p)
            close_drawing!(e)
        elseif near_first(e, p)                              # not the first vertex again
            e.status = "A zone needs at least 3 vertices."
        else
            q, kind = snap(e, p)
            push!(d, q); kind == :edge && push!(e.onedge, q)
        end
        return e
    end
    v = near_vertex(e, p)
    if isnothing(v)                                          # on an edge: insert a vertex there, and drag it
        snapshot = deepcopy(e.regions)
        v = insert_on_edge!(e.regions, p, tol(e))
        isnothing(v) && return e
        push!(e.history, snapshot); empty!(e.future)
    else
        checkpoint!(e)
    end
    e.drag = shared(e.regions, v)
    return e
end
function move!(e::Editor, x, z)
    p = [Float64(x), Float64(z)]
    if e.step == 1
        e.cursor = nothing
        if e.resizing
            s = nice_step(max(p[1], p[2], 1e-6))
            set_domain!(e, [max(s, round(p[1]/s)*s), max(s, round(p[2]/s)*s)])
        end
        e.readout = "x = $(round(p[1]; digits=3)) m, z = $(round(p[2]; digits=3)) m"
        return e
    end
    e.step == 2 || (e.cursor = nothing; e.readout = ""; return e)
    q, kind = closing(e, p) ? (copy(e.drawing[1]), :close) : snap(e, p; skip=e.drag)
    e.cursor = (q, kind)
    for (k, i) ∈ e.drag
        e.regions[k].vertices[i] = copy(q)
    end
    e.readout = "x = $(round(q[1]; digits=3)) m, z = $(round(q[2]; digits=3)) m" *
                (kind == :close ? " (click to close the zone)" : kind == :vertex ? " (on a vertex)" : kind == :edge ? " (on an edge)" : "")
    return e
end
function release!(e::Editor)
    e.resizing && (e.resizing = false; fit!(e))
    empty!(e.drag)
    return e
end
leave!(e::Editor) = (e.cursor = nothing; e.readout = ""; e)
function key!(e::Editor, key::AbstractString; cmd::Bool=false, shift::Bool=false)
    e.step == 2 || return e
    if cmd && (key == "y" || (key == "z" && shift))
        return redo!(e)
    elseif cmd && key == "z"
        isnothing(e.drawing) && return undo!(e)
        isempty(e.drawing) || (q = pop!(e.drawing); filter!(!=(q), e.onedge))
        return e
    end
    isnothing(e.drawing) && return e
    if key == "enter"
        close_drawing!(e)
    elseif key == "backspace"
        isempty(e.drawing) || (q = pop!(e.drawing); filter!(!=(q), e.onedge))
    elseif key == "escape"
        e.drawing = nothing; e.status = ""
    end
    return e
end
# the drawing as a zone: its edge-snapped vertices are inserted into the edges they lie on (shared vertices)
function close_drawing!(e::Editor)
    d = e.drawing
    length(d) < 3 && (e.status = "A zone needs at least 3 vertices."; return e)
    checkpoint!(e)
    foreach(q -> insert_on_edge!(e.regions, q, 1e-9*maximum(e.L); exact=true), e.onedge)
    empty!(e.onedge)
    push!(e.regions, DrawnRegion("zone $(length(e.regions) + 1)", d, Dict{String,Float64}()))
    e.drawing, e.status, e.selected = nothing, "", length(e.regions)
    return e
end

# ---- zones

new_zone!(e::Editor) = (e.drawing = Vector{Float64}[]; empty!(e.onedge); e.status = "Drawing: click vertices, close on the first one or press Enter."; e)
select_zone!(e::Editor, i) = (e.selected = clamp(Int(i), isempty(e.regions) ? 0 : 1, length(e.regions)); e)
function rename_zone!(e::Editor, name)
    i = e.selected; (i == 0 || name == e.regions[i].name) && return e
    checkpoint!(e); e.regions[i].name = String(name)
    return e
end
function set_prop!(e::Editor, key, text)
    i = e.selected; i == 0 && return e
    scale = PROPS[findfirst(t -> t[1] == key, PROPS)][3]
    v = tryparse(Float64, strip(text))
    get(e.regions[i].props, key, nothing) == (isnothing(v) ? nothing : v/scale) && return e
    checkpoint!(e)
    isnothing(v) ? delete!(e.regions[i].props, key) : (e.regions[i].props[key] = v/scale)
    return e
end
function swap_zone!(e::Editor, d)
    i, j = e.selected, e.selected + d
    (i == 0 || !(1 ≤ j ≤ length(e.regions))) && return e
    checkpoint!(e); e.regions[i], e.regions[j] = e.regions[j], e.regions[i]; e.selected = j
    return e
end
function delete_zone!(e::Editor)
    i = e.selected; i == 0 && return e
    checkpoint!(e); deleteat!(e.regions, i); select_zone!(e, i)
    return e
end
function restore!(e::Editor, from, to)
    isempty(from) && return e
    push!(to, deepcopy(e.regions)); s = pop!(from)
    empty!(e.regions); append!(e.regions, s); select_zone!(e, e.selected)
    return e
end
undo!(e::Editor) = restore!(e, e.history, e.future)
redo!(e::Editor) = restore!(e, e.future, e.history)

# ---- solver options: rows by key ("basis.which", "output.what.point.solid")

const PAGE_ORDER = ["basis", "stab", "material", "grf", "nonloc", "plot", "output", "output · point", "output · mesh", "dtype", "backend"]
option_sections() = sort(unique([option_page(r[1], r[2]) for r ∈ [option_rows(); option_lists()]]);
                         by=s -> something(findfirst(==(s), PAGE_ORDER), length(PAGE_ORDER) + 1))
option_key(sec, path) = join([sec; path], ".")
function option_entry(key)
    for r ∈ [option_rows(); option_lists()]
        option_key(r[1], r[2]) == key && return r
    end
    throw(ArgumentError("no solver option $key"))
end
# a list's values: plot.what offers the material-point variables only (the mesh ones need a `what` field)
list_values(sec, vals) = sec == "plot" ? filter(v -> haskey(v, :mpts), vals) : vals

# pick value j of a drop-down row (ticks it)
set_option!(e::Editor, key, j) = (r = option_entry(key); setoption!(e.options, r[1], r[2], r[3][Int(j)]); e)
# type a value into a numeric row (ticks it); an invalid text is ignored
function type_option!(e::Editor, key, text)
    r = option_entry(key)
    v = try option_value([r[1]; r[2]], text) catch; return e end
    setoption!(e.options, r[1], r[2], v)
    return e
end
# tick or untick a row (unticked: the default again)
function tick_option!(e::Editor, key, on::Bool)
    sec, path, _, def = option_entry(key)
    on ? setoption!(e.options, sec, path, something(getoption(e.options, sec, path), def)) : deleteoption!(e.options, sec, path)
    return e
end
# tick or untick value i of a list row (ticks the row)
function tick_value!(e::Editor, key, i, on::Bool)
    sec, path, vals, def = option_entry(key)
    vals = list_values(sec, vals)
    cur  = collect(Any, something(getoption(e.options, sec, path), def))   # Any: plot variables differ in type
    v    = vals[Int(i)]
    names = option_label.(cur)
    on ? (option_label(v) ∈ names || push!(cur, v)) : filter!(c -> option_label(c) != option_label(v), cur)
    setoption!(e.options, sec, path, [x for x ∈ vals if option_label(x) ∈ option_label.(cur)])   # in the listed order
    return e
end

# ---- run settings

function set_run!(e::Editor, key, text)
    t = strip(text)
    if key == "workflow"
        t ∈ WORKFLOWS && (e.run["workflow"] = String(t))
    elseif key == "fid"
        isempty(t) || (e.run["fid"] = String(t))
    elseif key == "lstps"
        n = tryparse(Int, t); !isnothing(n) && n ≥ 1 && (e.run["lstps"] = n)
    else
        x = tryparse(Float64, t); !isnothing(x) && x ≥ 0 && (e.run[key] = x)
    end
    return e
end

# ---- the case, and finishing

current_case(e::Editor) = Case(; L=copy(e.L), nel=copy(e.nel), ni=e.ni, width=e.width, nely=e.nely, bcs=copy(e.bcs),
                                material=e.base.material, thermal=e.base.thermal, solver=deepcopy(e.options), run=copy(e.run),
                                zones=[Zone(rg.name, Polygon(rg.vertices), merge(Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ rg.props), rg.extra))
                                       for rg ∈ e.regions])

# ---- cached: the preview points and the checks, recomputed only when what they depend on changes (a mouse
# move repaints the canvas, so a fine mesh would otherwise redo them at every move)

function cached!(f, e::Editor, name::Symbol, key::UInt)
    c = get(e.cache, name, nothing)
    !isnothing(c) && c[1] == key && return c[2]
    v = f()
    e.cache[name] = (key, v)
    return v
end
"""
    preview_points(editor) -> (x, z, zone)

The points `case_points` will create (`drawn_points`). While vertices are dragged, the last ones: they
follow on release.
"""
function preview_points(e::Editor)
    !isempty(e.drag) && haskey(e.cache, :points) && return e.cache[:points][2]
    return cached!(() -> drawn_points(e.L, e.nel, e.ni, e.regions), e, :points,
                   hash((e.L, e.nel, e.ni, [rg.vertices for rg ∈ e.regions])))
end
"`check_case` on the edited case (errors, warnings, info)."
function checks(e::Editor)
    key = hash((e.base, e.L, e.nel, e.ni, e.width, e.nely, e.bcs, e.options, e.run,
                [(rg.name, rg.vertices, rg.props, rg.extra) for rg ∈ e.regions]))
    return cached!(e, :checks, key) do
        try check_case(current_case(e)) catch err
            (; errors=[sprint(showerror, err)], warnings=String[], info=String[])
        end
    end
end
"""
    done!(editor) -> Bool

**Done**: the case, if it has no errors and (after a second Done) its warnings are accepted; whether
the window should close.
"""
function done!(e::Editor)
    c = checks(e)
    if !isempty(c.errors)
        e.status = join(c.errors, "\n")
    elseif !isempty(c.warnings) && !e.confirm
        e.status = join(c.warnings, "\n") * "\nClick Done again to continue anyway."; e.confirm = true
    else
        e.result, e.finished = current_case(e), true
        e.status = "Done: the case is back in Julia."
    end
    return e.finished
end
function cancel!(e::Editor)
    e.result, e.finished, e.status = nothing, true, "Cancelled."
    return e
end
function save!(e::Editor, file)
    try
        save_case(file, current_case(e)); e.status = "Saved $file"
    catch err
        e.status = "Not saved: " * sprint(showerror, err)
    end
    return e
end
function load!(e::Editor, file)
    endswith(lowercase(file), ".toml") || (e.status = "Not a .toml file: $file"; return e)
    c = try load_case(file) catch err
        e.status = "Not loaded: " * sprint(showerror, err); return e
    end
    try
        n = Editor(c; step=e.step)
        for f ∈ fieldnames(Editor)
            setfield!(e, f, getfield(n, f))
        end
        e.status = "Loaded $file"
    catch err
        e.status = "Not loaded: " * sprint(showerror, err)
    end
    return e
end
