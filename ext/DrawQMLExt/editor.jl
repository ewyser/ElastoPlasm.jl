# The drawing window as plain Julia state: what is drawn, and what each user action does. A GUI
# (`DrawQMLExt`) only shows it and forwards events, so everything here is testable without one.
# Coordinates are data coordinates (x, z in m) unless named `px`.

# a picture to trace the zones over (a scanned section, a map), drawn under the grid: its file, its size in pixels
# (reported by Qt once loaded; (0, 0) until then), its lower-left corner and width in metres (the height follows
# from the aspect), and its opacity. A drawing aid only: not part of the case.
mutable struct Backdrop
    file    :: String
    npx     :: Tuple{Int,Int}
    x0      :: Float64
    z0      :: Float64
    width   :: Float64
    opacity :: Float64
end
height(b::Backdrop) = b.npx[1] == 0 ? 0.0 : b.width*b.npx[2]/b.npx[1]

"""
    Editor(case; step=2)

The state of a drawing window editing `case`: the zones (polygons and circles), the domain, the boundary
conditions, the solver overrides and run settings, undo/redo snapshots, the polygon being drawn,
the vertices being dragged, the current step (1 Domain, 2 Zones, 3 Options, 4 Run), the options page,
the selected zone, the snapping toggles, the visible area (`view`, and the point grabbed while panning it), the
cursor, its coordinates (`readout`, drawn on the canvas) and the status line. `cache` keeps the preview points and the checks
until what they depend on changes (`preview_points`, `checks`).
`result` is the `Case` after **Create…** (`create!`), or `nothing` when cancelled; `finished` says
whether the window should close. The picture to trace over (`backdrop`, shown while `picture`) is kept when a
case is loaded.
"""
Base.@kwdef mutable struct Editor
    # from the case
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
    step     :: Int
    selected :: Int
    # the window's own state
    history  :: Vector{Vector{DrawnRegion}}                    = Vector{DrawnRegion}[]
    future   :: Vector{Vector{DrawnRegion}}                    = Vector{DrawnRegion}[]
    drawing  :: Union{Nothing,Vector{Vector{Float64}}}         = nothing
    onedge   :: Vector{Vector{Float64}}                        = Vector{Float64}[]   # vertices of the drawing snapped onto edges
    drag     :: Vector{Tuple{Int,Int}}                         = Tuple{Int,Int}[]    # (zone, vertex) moving together
    resizing :: Bool                                           = false   # the domain corner is being dragged
    grab     :: Union{Nothing,Vector{Float64}}                 = nothing   # panning: the point grabbed, kept under the cursor
    clicked  :: Union{Nothing,Tuple{Int,NTuple{4,Float64}}}    = nothing   # the zone pressed in, and the view then (unless panned: selected)
    section  :: String                                         = first(SECTIONS)
    grid     :: Bool                                           = true
    shapes   :: Bool                                           = true
    points   :: Bool                                           = true
    view     :: NTuple{4,Float64}                              = (0.0, 1.0, 0.0, 1.0)   # x0, x1, z0, z1
    cursor   :: Union{Nothing,Tuple{Vector{Float64},Symbol}}   = nothing
    status   :: String                                         = ""
    confirm  :: Bool                                           = false   # warnings shown: the next Create… goes through
    result   :: Union{Nothing,Case}                            = nothing
    finished :: Bool                                           = false
    readout  :: String                                         = ""      # the cursor's coordinates
    cache    :: Dict{Symbol,Tuple{UInt,Any}}                   = Dict{Symbol,Tuple{UInt,Any}}()   # name => (key, value)
    shape    :: Symbol                                         = :polygon   # the zone being drawn: :polygon, :rectangle, :triangle, :circle
    ring     :: Union{Nothing,Tuple{Int,Symbol}}               = nothing   # a circle being dragged: (zone, :centre or :radius)
    picked   :: Union{Nothing,Vector{Float64}}                 = nothing   # the vertex last clicked or dragged (its position)
    backdrop :: Union{Nothing,Backdrop}                        = nothing   # the picture under the grid …
    picture  :: Bool                                           = false     # … shown: the Domain step's choice, draw on a picture
    calib    :: Union{Nothing,Vector{Vector{Float64}}}         = nothing   # measuring a scale bar on it: its ends clicked (at most 2)
end
function Editor(case::Case; step::Integer=2)
    length(case.L) == 2 || throw(ArgumentError("the window draws a 2D (x–z) case; a width (domain step) extrudes it along y"))
    e = Editor(; base=case, regions=drawn_regions(case), L=copy(case.L), nel=copy(case.nel), h=minimum(case.L ./ case.nel),
               case.ni, case.width, case.nely, bcs=copy(case.bcs), options=config_dict(case.solver),
               run=Dict{String,Any}(string(k) => v for (k, v) ∈ pairs(case.run)), step=Int(step),
               selected=isempty(case.zones) ? 0 : 1)
    return fit!(e)
end

lattice(e::Editor) = (e.L[1]/(e.nel[1]*e.ni), e.L[2]/(e.nel[2]*e.ni))
tol(e::Editor)     = 0.015*max(e.view[2] - e.view[1], e.view[4] - e.view[3])   # pick and snap radius [m]: follows the zoom
"Zoom by `f` (< 1: in) about the point `(x, z)`, which stays where it is on screen."
function zoom!(e::Editor, x, z, f)
    x0, x1, z0, z1 = e.view
    e.view = (x + (x0 - x)*f, x + (x1 - x)*f, z + (z0 - z)*f, z + (z1 - z)*f)
    return e
end
"Move the visible area by `d` (x, z) [m]."
function pan!(e::Editor, d)
    x0, x1, z0, z1 = e.view
    e.view = (x0 + d[1], x1 + d[1], z0 + d[2], z1 + d[2])
    return e
end
"""
    wheel!(editor, x, z, angle, shift; zoom=false)

The wheel at `(x, z)`: a trackpad's two-finger scroll (`shift`, how far the content moves, x and z [m], non-zero)
pans the view; a mouse wheel (`angle`, in eighths of a degree, 120 a notch), or any scroll with `zoom` (Cmd/Ctrl
held), zooms about the pointer.
"""
function wheel!(e::Editor, x, z, angle, shift; zoom::Bool=false)
    !zoom && any(!=(0), shift) && return pan!(e, .-shift)
    return zoom!(e, x, z, 0.9^(angle/120))
end
"Fit the visible area to the domain, with a small margin."
function fit!(e::Editor)
    e.view = (-0.02e.L[1], 1.02e.L[1], -0.04e.L[2], 1.04e.L[2])
    return e
end
# an undo step before an edit, or with the zones as they were before it (`snapshot`, when the edit may not happen)
function checkpoint!(e::Editor, snapshot=deepcopy(e.regions))
    push!(e.history, snapshot); length(e.history) > 200 && popfirst!(e.history); empty!(e.future)
    return e
end
outside(e::Editor) = any(v -> !(0 ≤ v[1] ≤ e.L[1] && 0 ≤ v[2] ≤ e.L[2]), (v for rg ∈ e.regions for v ∈ rg.vertices)) ||
                     any(rg -> iscircle(rg) && !all(0 .≤ rg.centre .- rg.radius .&& rg.centre .+ rg.radius .≤ e.L), e.regions)
outside!(e::Editor) = (e.status = outside(e) ? "Some zones lie outside the domain." : ""; e)
# a circle's handles: its centre, and the point east of it that sets the radius
radius_handle(rg::DrawnRegion) = [rg.centre[1] + rg.radius, rg.centre[2]]
function near_handle(e::Editor, p)                           # (zone, :centre or :radius) of a circle under p, or nothing
    for k ∈ reverse(eachindex(e.regions))
        rg = e.regions[k]; iscircle(rg) || continue
        dist(p, rg.centre) < tol(e) && return (k, :centre)
        dist(p, radius_handle(rg)) < tol(e) && return (k, :radius)
    end
    return nothing
end
snap(e::Editor, p; skip=()) = snap_point(p, e.regions, lattice(e), e.L, tol(e); grid=e.grid, shapes=e.shapes, skip)
# the zone being drawn: is `p` on its first vertex (a click there closes it, from 3 vertices on)?
near_first(e::Editor, p) = !isnothing(e.drawing) && !isempty(e.drawing) && dist(p, e.drawing[1]) < tol(e)
closing(e::Editor, p) = e.shape == :polygon && near_first(e, p) && length(e.drawing) ≥ 3
near_vertex(e::Editor, p) = (for rg ∈ reverse(e.regions), v ∈ rg.vertices
                                 dist(p, v) < tol(e) && return v
                             end; nothing)

# ---- the domain

function set_domain!(e::Editor, L; h=e.h)
    e.h, e.L, e.nel = h, Float64.(collect(L)), nel_for(L, h)
    isnothing(e.width) || (e.nely = nel_for([e.width], h)[1])
    return outside!(e)
end
function set_width!(e::Editor, text)
    x = number(text)
    e.width = isnothing(x) || x ≤ 0 ? nothing : x
    e.nely  = isnothing(e.width) ? nothing : nel_for([e.width], e.h)[1]
    return e
end
"Set the domain field `name` (`Lx`, `Lz`, the element size `h`, or the `width` along y) to the number in `text`."
function set_domain_field!(e::Editor, name, text)
    name == "width" && return set_domain!(set_width!(e, text), e.L)
    x = number(text); (isnothing(x) || x ≤ 0) && return e
    name == "h" ? set_domain!(e, e.L; h=x) : set_domain!(e, name == "Lx" ? [x, e.L[2]] : [e.L[1], x])
    return fit!(e)
end
set_bc!(e::Editor, side, kind) = (s = findfirst(t -> t[1] == side, SIDES); e.bcs[SIDES[s][2], SIDES[s][3]] = Symbol(kind); e)
function set_step!(e::Editor, k)
    e.step, e.drawing, e.confirm, e.cursor = Int(k), nothing, false, nothing
    fit!(e)
    return e
end

# ---- mouse and keys (data coordinates)

function press!(e::Editor, x, z; right::Bool=false, pan::Bool=false)
    p = [Float64(x), Float64(z)]
    pan && (e.grab = p; return e)                           # the middle button: pan, whatever is under it
    if !isnothing(e.calib) && length(e.calib) < 2 && !right  # a scale bar's end on the picture (not snapped)
        push!(e.calib, p)
        e.status = length(e.calib) == 1 ? "Scale bar: click its other end." : "Scale bar: type its length in metres, then Enter."
        return e
    end
    if e.step == 1                                           # the domain: grab the corner handle, or pan
        right && return e
        dist(p, e.L) < 2tol(e) ? (e.resizing = true) : (e.grab = p)
        return e
    end
    e.step == 2 || (right || (e.grab = p); return e)          # options, run: the canvas only pans
    if right                                                 # delete a vertex
        v = isnothing(e.drawing) ? near_vertex(e, p) : nothing
        isnothing(v) || remove_at!(e, v)
        return e
    end
    if !isnothing(e.drawing) && e.shape == :circle          # drawing a circle: its centre, then a point on it
        q, _ = snap(e, p)
        if isempty(e.drawing)
            push!(e.drawing, q)
        else
            r = dist(q, e.drawing[1])
            r > 0 || return e
            checkpoint!(e)
            push!(e.regions, circle_region("zone $(length(e.regions) + 1)", e.drawing[1], r, Dict{String,Float64}()))
            e.drawing, e.status, e.selected = nothing, "", length(e.regions)
            cut_earlier!(e)
        end
        return e
    end
    if !isnothing(e.drawing) && e.shape == :rectangle       # drawing a rectangle: a corner, then the opposite one
        q, _ = snap(e, p)
        if isempty(e.drawing)
            push!(e.drawing, q)
        else
            a = e.drawing[1]
            (q[1] == a[1] || q[2] == a[2]) && (e.status = "A rectangle needs two opposite corners."; return e)
            e.drawing = [a, [q[1], a[2]], q, [a[1], q[2]]]
            append!(e.onedge, e.drawing)                     # every corner on another zone's edge splits it (shared)
            close_drawing!(e)
        end
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
            e.shape == :triangle && length(d) == 3 && close_drawing!(e)   # a triangle closes on its third corner
        end
        return e
    end
    h = near_handle(e, p)
    isnothing(h) || (checkpoint!(e); e.ring = h; return e)  # a circle's centre or radius handle: drag it
    v = near_vertex(e, p)
    if isnothing(v)                                          # on an edge: insert a vertex there, and drag it
        snapshot = deepcopy(e.regions)
        v = insert_on_edge!(e.regions, p, tol(e))
        if isnothing(v)                                      # on no vertex or edge: pan, or select the zone clicked in
            e.grab = p
            k = findlast(rg -> inside(shape_of(rg), (p[1], p[2])), e.regions)   # the topmost, whose constants a point there gets
            isnothing(k) || (e.clicked = (k, e.view))
            return e
        end
        checkpoint!(e, snapshot)
    else
        checkpoint!(e)
    end
    e.drag = shared(e.regions, v)
    e.picked = copy(v); e.selected = last(e.drag)[1]               # its row in the panel; its (topmost) zone selected
    return e
end
function move!(e::Editor, x, z)
    p = [Float64(x), Float64(z)]
    isnothing(e.grab) || return pan!(e, e.grab .- p)         # panning: the grabbed point follows the cursor
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
    isempty(e.drag) || (e.picked = copy(q))
    if !isnothing(e.ring)                                    # a circle: move its centre, or set its radius
        k, what = e.ring; rg = e.regions[k]
        what == :centre ? (rg.centre = copy(q)) :
                          (rg.radius = max(dist(q, rg.centre), minimum(lattice(e))/2))
    end
    e.readout = "x = $(round(q[1]; digits=3)) m, z = $(round(q[2]; digits=3)) m" *
                (kind == :close ? " (click to close the zone)" : kind == :vertex ? " (on a vertex)" : kind == :edge ? " (on an edge)" : "")
    return e
end
function release!(e::Editor)
    e.resizing && (e.resizing = false; fit!(e))
    if !isnothing(e.clicked)                                 # a click inside a zone, not a pan (within a few pixels), selects it
        k, view = e.clicked
        maximum(abs.(e.view .- view)) < 0.2tol(e) && (e.selected = k)
    end
    empty!(e.drag); e.ring, e.grab, e.clicked = nothing, nothing, nothing
    return e
end
leave!(e::Editor) = (e.cursor = nothing; e.readout = ""; e)
function key!(e::Editor, key::AbstractString; cmd::Bool=false, shift::Bool=false)
    key == "escape" && !isnothing(e.calib) && (e.calib = nothing; e.status = ""; return e)
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
        e.shape == :polygon && close_drawing!(e)
    elseif key == "backspace"
        isempty(e.drawing) || (q = pop!(e.drawing); filter!(!=(q), e.onedge))
    elseif key == "escape"
        e.drawing = nothing; e.status = ""
    end
    return e
end
# the drawing as a zone: its edge-snapped vertices are inserted into the edges they lie on (shared vertices), and it
# is cut out of the earlier zones
function close_drawing!(e::Editor)
    d = e.drawing
    length(d) < 3 && (e.status = "A zone needs at least 3 vertices."; return e)
    checkpoint!(e)
    foreach(q -> insert_on_edge!(e.regions, q, 1e-9*maximum(e.L); exact=true), e.onedge)
    empty!(e.onedge)
    push!(e.regions, DrawnRegion("zone $(length(e.regions) + 1)", d, Dict{String,Float64}()))
    e.drawing, e.status, e.selected = nothing, "", length(e.regions)
    return cut_earlier!(e)
end
"""
    cut_earlier!(editor)

Cut the zone just drawn (the last) out of the earlier zones it overlaps (`cut`): each zone's outline is then where
its points are (a point belongs to the last zone containing it, so the points do not change, but for a cut circle's
next to its arc). An earlier zone may be left in several parts or with holes (several outlines); one covered
entirely is removed; a circle that is cut becomes a polygon (`zone_rings`). Where the new zone's outline crosses
theirs it gets a vertex too, so that the boundary they now share moves together. Part of the same undo step as
drawing the zone.
"""
function cut_earlier!(e::Editor)
    new, scale = last(e.regions), maximum(e.L)
    cuts, gone = String[], Int[]
    for k ∈ 1:length(e.regions)-1
        rg = e.regions[k]
        overlap(bbox(rg), bbox(new)) || continue
        rs = cut(zone_rings(rg), [zone_rings(new)], scale)
        isnothing(rs) && continue
        isempty(rs) && (push!(gone, k); continue)
        rg.vertices, rg.sizes, rg.centre, rg.radius = reduce(vcat, rs), length.(rs), nothing, 0.0
        push!(cuts, rg.name)
        iscircle(new) && continue
        for v ∈ rg.vertices                                  # where the outlines cross: a vertex of the new zone too
            v ∈ new.vertices || insert_on_edge!([new], v, 1e-9scale; exact=true)
        end
    end
    msg = String[]
    isempty(cuts) || push!(msg, "Cut out of $(join(cuts, ", ")).")
    isempty(gone) || push!(msg, "Covered by it, removed: $(join((e.regions[k].name for k ∈ gone), ", ")).")
    deleteat!(e.regions, gone)
    e.status, e.selected = join(msg, " "), length(e.regions)
    return e
end

# ---- zones

const NEW_ZONE = Dict(:polygon   => "Draw: click vertices, close on the first one or press Enter.",
                      :rectangle => "Rectangle: click a corner, then the opposite one.",
                      :triangle  => "Triangle: click its three corners.",
                      :circle    => "Circle: click its centre, then a point on it.")
"Start a new zone of `shape`: `:polygon` (drawn vertex by vertex), `:rectangle`, `:triangle` or `:circle`."
function new_zone!(e::Editor, shape::Symbol=:polygon)
    haskey(NEW_ZONE, shape) || throw(ArgumentError("no zone shape $shape; one of $(collect(keys(NEW_ZONE)))"))
    e.drawing, e.shape, e.status = Vector{Float64}[], shape, NEW_ZONE[shape]; empty!(e.onedge)
    return e
end
"Set the centre's `x` or `z`, or the radius `r`, of the selected circle to the number in `text` (otherwise ignored)."
function set_circle!(e::Editor, key, text)
    k = e.selected; (k == 0 || !iscircle(e.regions[k])) && return e
    x = number(text); (isnothing(x) || (key == "r" && x ≤ 0)) && return e
    rg = e.regions[k]
    (key == "r" ? rg.radius : rg.centre[key == "x" ? 1 : 2]) == x && return e
    checkpoint!(e)
    key == "r" ? (rg.radius = x) : (rg.centre[key == "x" ? 1 : 2] = x)
    return outside!(e)
end
select_zone!(e::Editor, i) = (e.selected = clamp(Int(i), isempty(e.regions) ? 0 : 1, length(e.regions)); e)
function rename_zone!(e::Editor, name)
    i = e.selected; (i == 0 || name == e.regions[i].name) && return e
    checkpoint!(e); e.regions[i].name = String(name)
    return e
end
function set_prop!(e::Editor, key, text)
    i = e.selected; i == 0 && return e
    scale = PROPS[findfirst(t -> t[1] == key, PROPS)][3]
    v = number(text)
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

# ---- the picture to trace over. The Domain step starts with a choice: a domain sized by hand, or a picture to draw
# on (`use_picture!`), which has its own panel: load the picture, give its scale (a scale bar: two clicks and its
# length, or its width) and place it; the domain follows it, from the origin to its far corner (`fit_domain!`).
# Its lower-left corner (x0, z0) is where it sits in the domain (0, 0 by default; move it to crop the picture).

"Draw on the picture (`on`: shown, its panel in the Domain step), or on a domain sized by hand (the picture kept, hidden)."
function use_picture!(e::Editor, on::Bool)
    e.picture, e.calib = on, nothing
    e.status = on && isnothing(e.backdrop) ? "Choose a picture to draw on: a scanned section, a map." : ""
    return e
end
"Show the picture `file` under the grid, across the domain's width until its scale is given."
function set_backdrop!(e::Editor, file)
    e.backdrop = Backdrop(String(file), (0, 0), 0.0, 0.0, e.L[1], 0.6)
    e.picture, e.calib = true, nothing
    e.status = "Loading the picture…"
    return e
end
function remove_backdrop!(e::Editor)
    e.backdrop, e.calib = nothing, nothing
    return e
end
"The picture's size in pixels, reported by Qt once it has loaded it: the domain is fitted to the picture."
function backdrop_size!(e::Editor, w, h)
    isnothing(e.backdrop) && return e
    e.backdrop.npx = (Int(w), Int(h))
    fit_domain!(e)
    e.status = "Picture loaded: measure a scale bar on it, or type its width."   # after the fit, which resets it
    return e
end
"""
    set_backdrop_field!(editor, key, text)

Set the picture's `width` [m] (its scale), `x0`/`z0` (its lower-left corner in the domain) or `opacity` (0–1) to
the number in `text`; other text is ignored. The domain is fitted to the picture again, unless the opacity changed.
"""
function set_backdrop_field!(e::Editor, key, text)
    b = e.backdrop; isnothing(b) && return e
    x = number(text); isnothing(x) && return e
    if key == "width" && x > 0
        b.width = x
    elseif key == "x0"
        b.x0 = x
    elseif key == "z0"
        b.z0 = x
    elseif key == "opacity"
        b.opacity = clamp(x, 0.0, 1.0)
        return e
    else
        return e
    end
    return fit_domain!(e)
end
"""
    fit_domain!(editor)

Fit the domain to the picture, from the origin to its far corner, in whole square elements of the current size
(the nearest number of them): a measured scale gives an arbitrary picture size, which would make the cells
slightly elongated.
"""
function fit_domain!(e::Editor)
    b = e.backdrop
    (isnothing(b) || b.npx[1] == 0) && return e
    L = [b.x0 + b.width, b.z0 + height(b)]
    all(>(0), L) || return e
    set_domain!(e, e.h .* max.(1, round.(L ./ e.h))); fit!(e)
    return e
end
"Measure a scale bar on the picture: the next two clicks are its ends, then its length is typed (`apply_scale!`)."
function scale_bar!(e::Editor)
    isnothing(e.backdrop) && return e
    e.calib = Vector{Float64}[]
    e.status = "Scale bar: click one end of it on the picture."
    return e
end
"""
    apply_scale!(editor, text)

Scale the picture so that the scale bar clicked on it is as long as the number in `text` [m] (other text is
ignored): a uniform scale about the picture's lower-left corner (the aspect is kept), from the ratio of that
length to the clicked one; the domain is then fitted to the picture.
"""
function apply_scale!(e::Editor, text)
    b, c, len = e.backdrop, e.calib, number(text)
    (isnothing(b) || isnothing(c) || length(c) < 2 || isnothing(len) || len ≤ 0) && return e
    d = dist(c[2], c[1])
    d > 0 || (e.status = "The two ends must differ."; return e)
    b.width *= len/d
    e.calib = nothing
    fit_domain!(e); e.status = "Picture scaled, domain fitted to it."
    return e
end

# ---- the vertex last clicked or dragged (`picked`), typed exactly in the Zones step. A vertex shared with other
# zones moves, or is deleted, in all of them, as when dragging or right-clicking it on the canvas.

# the picked vertex, or nothing once no zone has it (undo, deletion)
picked_vertex(e::Editor) = isnothing(e.picked) || isempty(shared(e.regions, e.picked)) ? nothing : e.picked
"Set coordinate `axis` (\"x\" or \"z\") of the picked vertex to the number in `text` (otherwise ignored)."
function set_picked!(e::Editor, axis, text)
    v = e.picked; isnothing(v) && return e
    at = shared(e.regions, v); isempty(at) && (e.picked = nothing; return e)   # gone (undo, deletion)
    x = number(text); isnothing(x) && return e
    w = copy(v); w[axis == "x" ? 1 : 2] = x
    w == v && return e
    checkpoint!(e)
    foreach(((k, i),) -> e.regions[k].vertices[i] = copy(w), at)
    e.picked = w
    return outside!(e)
end
"Delete the picked vertex (a zone keeps at least 3)."
remove_picked!(e::Editor) = (isnothing(e.picked) || (remove_at!(e, e.picked); isempty(shared(e.regions, e.picked)) && (e.picked = nothing)); e)
function remove_at!(e::Editor, v)                            # the vertex at `v`, from every zone that has it
    snapshot = deepcopy(e.regions)
    delete_vertex!(e.regions, v) ? checkpoint!(e, snapshot) : (e.status = "A zone keeps at least 3 vertices.")
    return e
end

# ---- solver options: rows by key ("basis.which", "output.what.point.solid")

const PAGE_ORDER = ["basis", "stab", "material", "grf", "nonloc", "plot", "output", "dtype", "backend"]
const SECTIONS   = sort(unique(first.(OPTIONS));                                  # the pages of the section menu, one per section
                        by=s -> something(findfirst(==(s), PAGE_ORDER), length(PAGE_ORDER) + 1))
option_key(sec, path) = join([sec; path], ".")
function option_entry(key)
    i = findfirst(r -> option_key(r[1], r[2]) == key, OPTIONS)
    isnothing(i) && throw(ArgumentError("no solver option $key"))
    return OPTIONS[i]
end

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
    cur = something(getoption(e.options, sec, path), def)
    new = on ? union(cur, [vals[Int(i)]]) : setdiff(cur, [vals[Int(i)]])
    setoption!(e.options, sec, path, [v for v ∈ vals if v ∈ new])                # in the listed order
    return e
end

# ---- run settings

function set_run!(e::Editor, key, text)
    t = strip(text)
    if key == "workflow"
        t ∈ WORKFLOWS && (e.run["workflow"] = String(t))
    elseif key == "lstps"
        n = tryparse(Int, t); !isnothing(n) && n ≥ 1 && (e.run["lstps"] = n)
    else
        x = number(t); !isnothing(x) && x ≥ 0 && (e.run[key] = x)
    end
    return e
end

# ---- the case, and finishing

current_case(e::Editor) = Case(; L=copy(e.L), nel=copy(e.nel), ni=e.ni, width=e.width, nely=e.nely, bcs=copy(e.bcs),
                                material=e.base.material, thermal=e.base.thermal, solver=deepcopy(e.options), run=copy(e.run),
                                zones=[Zone(rg.name, shape_of(rg), merge(Dict{Symbol,Any}(Symbol(k) => v for (k, v) ∈ rg.props), rg.extra))
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
    (!isempty(e.drag) || !isnothing(e.ring)) && haskey(e.cache, :points) && return e.cache[:points][2]
    return cached!(() -> drawn_points(e.L, e.nel, e.ni, e.regions), e, :points,
                   hash((e.L, e.nel, e.ni, [(rg.vertices, rg.sizes, rg.centre, rg.radius) for rg ∈ e.regions])))
end
"`check_case` on the edited case (errors, warnings, info)."
function checks(e::Editor)
    key = hash((e.base, e.L, e.nel, e.ni, e.width, e.nely, e.bcs, e.options, e.run,
                [(rg.name, rg.vertices, rg.sizes, rg.centre, rg.radius, rg.props, rg.extra) for rg ∈ e.regions]))
    return cached!(e, :checks, key) do
        try check_case(current_case(e)) catch err
            (; errors=[sprint(showerror, err)], warnings=String[], info=String[])
        end
    end
end
"""
    ready!(editor) -> Bool

**Create…**, before its folder dialog: whether the case may be created, i.e. it has no errors and its
warnings, shown on a first click, are accepted by a second.
"""
function ready!(e::Editor)
    c = checks(e)
    if !isempty(c.errors)
        e.status = join(c.errors, "\n"); return false
    elseif !isempty(c.warnings) && !e.confirm
        e.status = join(c.warnings, "\n") * "\nClick Create… again to continue anyway."; e.confirm = true; return false
    end
    return true
end
"""
    create!(editor, folder) -> Bool

**Create…** in `folder` (made if missing): the case's runs go there (`run.folder`),
the case is written there as `<name>_case.toml`, and the window closes with the case as its result.
"""
function create!(e::Editor, folder)
    dir  = rstrip(normpath(abspath(String(folder))), ['/', '\\'])
    name = basename(dir)
    isempty(name) && (e.status = "Choose a folder for the case."; return false)
    e.run["folder"] = dir
    c = current_case(e)
    try
        mkpath(dir); save_case(joinpath(dir, "$(name)_case.toml"), c)
    catch err
        e.status = "Not created: " * sprint(showerror, err); return false
    end
    e.result, e.finished, e.status = c, true, "Created in $dir."
    return true
end
function cancel!(e::Editor)                                  # once: closing the window after Create… cancels nothing
    e.finished || ((e.result, e.finished, e.status) = (nothing, true, "Cancelled."))
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
            f ∈ (:backdrop, :picture) || setfield!(e, f, getfield(n, f))   # the picture stays: the case is drawn on it
        end
        e.status = "Loaded $file"
    catch err
        e.status = "Not loaded: " * sprint(showerror, err)
    end
    return e
end
