# The parts of the drawing window that do not depend on the GUI: zone constants, points preview,
# snapping and edge editing, option rows. Used by the `Editor` (editor.jl).

const COLORS = ["#2f6fde", "#e0782b", "#2a9d6f", "#c43c6a", "#7a5cc7", "#b58b14", "#1f8fa6", "#6b7280"]
# the zone constants the window edits: key, label, display unit per SI unit (angles stored in radians)
const PROPS  = [("c0", "c₀ [kPa]", 1e-3), ("cr", "cᵣ [kPa]", 1e-3), ("phi", "φ [°]", 180/π), ("T", "T [K]", 1.0),
                ("E", "E [MPa]", 1e-6), ("nu", "ν [–]", 1.0), ("rho0", "ρ₀ [kg/m³]", 1.0)]
const GUI_KEYS = first.(PROPS)
const SIDES  = [("left", 1, 1), ("right", 1, 2), ("bottom", 2, 1), ("top", 2, 2)]   # (name, axis, lower/upper) of bcs

"""
    drawn_points(L, nel, ni, regions) -> (x, z, region)

The 2D material points `case_points` creates for these zones (sub-cell centres, as `mpts_populate`)
and the index of the region each belongs to, by the zones' own `inside` rule; the window previews them.
"""
function drawn_points(L, nel, ni, regions)
    hx, hz = L[1]/(nel[1]*ni), L[2]/(nel[2]*ni)
    x = [(i - 0.5)*hx for i ∈ 1:nel[1]*ni, _ ∈ 1:nel[2]*ni]
    z = [(j - 0.5)*hz for _ ∈ 1:nel[1]*ni, j ∈ 1:nel[2]*ni]
    shapes = shape_of.(regions)
    r = [something(findlast(s -> inside(s, (x[k], z[k])), shapes), 0) for k ∈ eachindex(x)]
    return x[r .> 0], z[r .> 0], r[r .> 0]
end

# a zone as the window edits it: a polygon (its vertices) or a circle (a `Ball`: its centre and radius, and
# no vertices), the constants of `PROPS` (SI, by key), and the zone's other properties (velocity, Hp, n₀,
# cₚ, k) carried through unchanged
mutable struct DrawnRegion
    name     :: String
    vertices :: Vector{Vector{Float64}}
    props    :: Dict{String,Float64}
    extra    :: Dict{Symbol,Any}
    centre   :: Union{Nothing,Vector{Float64}}               # a circle: its centre (x, z) …
    radius   :: Float64                                      # … and radius
end
DrawnRegion(name, vertices, props, extra=Dict{Symbol,Any}()) = DrawnRegion(name, vertices, props, extra, nothing, 0.0)
circle_region(name, centre, radius, props, extra=Dict{Symbol,Any}()) =
    DrawnRegion(name, Vector{Float64}[], props, extra, Float64.(collect(centre)), Float64(radius))
iscircle(rg::DrawnRegion) = !isnothing(rg.centre)
shape_of(rg::DrawnRegion) = iscircle(rg) ? Ball(rg.centre, rg.radius) : Polygon(rg.vertices)

# the window's regions from a case's zones: polygons, (x, z) blocks as their rectangle, (x, z) balls as circles
function drawn_regions(case::Case)
    map(case.zones) do z
        s = z.shape
        props = Dict{String,Float64}(string(k) => x for (k, x) ∈ z.props if string(k) ∈ GUI_KEYS)
        extra = Dict{Symbol,Any}(k => x for (k, x) ∈ z.props if string(k) ∉ GUI_KEYS)
        if s isa Ball && length(s.centre) == 2
            return circle_region(z.name, s.centre, s.radius, props, extra)
        end
        vs = s isa Polygon ? deepcopy(s.vertices) :
             s isa Block && length(s.lo) == 2 ? [[s.lo[1], s.lo[2]], [s.hi[1], s.lo[2]], [s.hi[1], s.hi[2]], [s.lo[1], s.hi[2]]] :
             throw(ArgumentError("the window draws polygons, (x, z) blocks and (x, z) balls; zone \"$(z.name)\" is a $(nameof(typeof(s)))"))
        return DrawnRegion(z.name, vs, props, extra)
    end
end

# closest point to `p` on the segment `a`–`b`
function closest_on_segment(p, a, b)
    d = b .- a
    t = iszero(d[1]^2 + d[2]^2) ? 0.0 : clamp(((p[1] - a[1])*d[1] + (p[2] - a[2])*d[2])/(d[1]^2 + d[2]^2), 0.0, 1.0)
    return a .+ t .* d
end
edges(rg) = ((i, mod1(i + 1, length(rg.vertices))) for i ∈ eachindex(rg.vertices))

"""
    snap_point(p, regions, step, L, tol; grid=true, shapes=true, skip=()) -> (point, kind)

Where a click at `p` lands: on an existing vertex within `tol` (`kind = :vertex`), else on an existing
edge within `tol` (`:edge`), else on the sub-cell lattice `step` (`:grid`), else at `p` (`:free`);
clamped to the domain. `skip` lists `(region, vertex)` pairs left out, with their edges (the vertices
being dragged).
"""
function snap_point(p, regions, step, L, tol; grid::Bool=true, shapes::Bool=true, skip=())
    p = [Float64(p[1]), Float64(p[2])]
    clampL(q) = [min(max(q[1], 0.0), L[1]), min(max(q[2], 0.0), L[2])]   # max(-0.0, 0.0) == 0.0
    if shapes
        best, d = nothing, tol
        for (k, rg) ∈ enumerate(regions), (i, v) ∈ enumerate(rg.vertices)
            (k, i) ∈ skip && continue
            (h = hypot(p[1] - v[1], p[2] - v[2])) < d && ((best, d) = (v, h))
        end
        isnothing(best) || return copy(best), :vertex
        for (k, rg) ∈ enumerate(regions), (i, j) ∈ edges(rg)
            ((k, i) ∈ skip || (k, j) ∈ skip) && continue
            q = closest_on_segment(p, rg.vertices[i], rg.vertices[j])
            (h = hypot(p[1] - q[1], p[2] - q[2])) < d && ((best, d) = (q, h))
        end
        isnothing(best) || return clampL(best), :edge
    end
    grid && return clampL([round(p[1]/step[1])*step[1], round(p[2]/step[2])*step[2]]), :grid
    return clampL(p), :free
end

# all (region, vertex) pairs at the position `v`: a vertex shared by several regions
shared(regions, v) = [(k, i) for (k, rg) ∈ enumerate(regions) for (i, w) ∈ enumerate(rg.vertices) if w == v]

"""
    insert_on_edge!(regions, p, tol; exact=false) -> Union{Nothing,Vector{Float64}}

Insert a vertex on the edge nearest to `p` (within `tol`), at the closest point, in every region that
has this edge, so a shared boundary stays shared. Returns the new vertex, or `nothing`. With
`exact = true` the vertex is `p` itself (a vertex already snapped onto the edge, which a new zone
keeps: projecting it again would move it by round-off and break the shared boundary).
"""
function insert_on_edge!(regions, p, tol; exact=false)
    best, d = nothing, tol
    for rg ∈ regions, (i, j) ∈ edges(rg)
        a, b = rg.vertices[i], rg.vertices[j]
        q = closest_on_segment(p, a, b)
        (h = hypot(p[1] - q[1], p[2] - q[2])) < d && q != a && q != b && ((best, d) = ((q, a, b), h))
    end
    isnothing(best) && return nothing
    q, a, b = best
    exact && (q = p)
    for rg ∈ regions                                           # (two loops: `break` leaves only the inner one)
        for (i, j) ∈ collect(edges(rg))
            e = (rg.vertices[i], rg.vertices[j])
            (e == (a, b) || e == (b, a)) && (insert!(rg.vertices, i + 1, copy(q)); break)
        end
    end
    return q
end

"""
    delete_vertex!(regions, v) -> Bool

Remove the vertex at `v` from every region that has it and keeps at least 3 vertices.
"""
function delete_vertex!(regions, v)
    done = false
    for rg ∈ regions
        i = findfirst(==(v), rg.vertices)
        (isnothing(i) || length(rg.vertices) ≤ 3) && continue
        deleteat!(rg.vertices, i); done = true
    end
    return done
end

# a round step for dragging the domain size: a tenth of its order of magnitude
nice_step(x) = 10.0^floor(log10(max(x, 1e-9)))/10
# elements along each axis for the element size h (square cells as far as L allows)
nel_for(L, h) = max.(1, round.(Int, L ./ h))

"""
    option_rows() -> Vector{(section, path, values, default)}
    option_lists() -> Vector{(section, path, values, default)}

The solver options the window offers, from `get_option()` (labels and values) and `get_default()`
(defaults): `option_rows` the single-value fields, one row each, `option_lists` the multi-selections
(`plot.what`, the `output.what` lists), a checkbox per value. Left out: `dtype.T0`/`dtype.bits` (set
from `dtype.precision` by `get_solver`) and `solution` (set by the run's workflow,
`case_solver_options`). `path` is the field's key path below the section.
"""
function option_entries()
    rows, opt, def = [], get_option(), get_default()
    function walk(o, d, sec, path)
        if o isa Tuple && length(o) == 2 && o[1] isa AbstractString
            (sec == "dtype" && path ∈ (["T0"], ["bits"])) || sec == "solution" || push!(rows, (sec, path, o[2], d))
        elseif o isa NamedTuple
            for k ∈ keys(o)
                walk(getfield(o, k), getfield(d, k), sec, [path; string(k)])
            end
        end
    end
    for k ∈ keys(opt)
        walk(getfield(opt, k), getfield(def, k), string(k), String[])
    end
    return rows
end
option_rows()  = filter(r -> !multi_select([r[1]; r[2]]), option_entries())
option_lists() = filter(r -> multi_select([r[1]; r[2]]), option_entries())
# the page of the section menu an option is on: its section, or for nested lists (`output.what.point.…`)
# the section and the group, e.g. "output · point"
option_page(sec, path) = length(path) ≥ 3 ? "$sec · $(path[2])" : sec
# the override at (section, path) in the nested options Dict, or nothing; and setting/removing it
function getoption(options, sec, path)
    d = get(options, sec, nothing)
    isempty(path) && return d
    for k ∈ path[1:end-1]
        d isa AbstractDict || return nothing
        d = get(d, k, nothing)
    end
    return d isa AbstractDict ? get(d, path[end], nothing) : nothing
end
function setoption!(options, sec, path, v)
    isempty(path) && return (options[sec] = v)
    d = get!(options, sec, Dict{String,Any}())
    for k ∈ path[1:end-1]
        d = get!(d, k, Dict{String,Any}())
    end
    d[path[end]] = v
end
function deleteoption!(options, sec, path)
    isempty(path) && return delete!(options, sec)
    chain = [options]
    for k ∈ [sec; path[1:end-1]]
        haskey(chain[end], k) || return
        push!(chain, chain[end][k])
    end
    delete!(chain[end], path[end])
    for (parent, k) ∈ zip(reverse(chain[1:end-1]), reverse([sec; path[1:end-1]]))   # drop emptied tables
        isempty(parent[k]) ? delete!(parent, k) : break
    end
end

# a number or number list as text (and back, `parse_numbers`)
fmt(v) = v isa AbstractVector ? join(v, ", ") : string(v)
