# Cutting a zone by the zones drawn after it. The difference's boundary is traced from the arrangement of all their
# edges, split where they cross or touch: a piece of edge is kept when the difference lies on exactly one side of it,
# probed with the even–odd rule `case_points` uses. Probing sides, rather than following crossings, keeps shared edges
# and vertices (the window's snapped boundaries) exact, where clipping algorithms need special cases for them.

# the signed angle from direction u to v (left turns positive)
function turn(u, v)
    return atan(u[1]*v[2] - u[2]*v[1], u[1]*v[1] + u[2]*v[2])
end
# the index of point p in `pts`, merged with one closer than `tol` (added if there is none)
function point_index!(pts, p, tol)
    k = findfirst(q -> dist(p, q) < tol, pts)
    isnothing(k) || return k
    push!(pts, p)
    return length(pts)
end

# where segment c–d meets segment a–b, as parameters along a–b (0 at a, 1 at b): a crossing or a touching end, or
# the ends of a collinear overlap, within `tol` [m]
function meets(a, b, c, d, tol)
    r, s, w = b .- a, d .- c, c .- a
    lr, ls = hypot(r...), hypot(s...)
    den = r[1]*s[2] - r[2]*s[1]
    ts = Float64[]
    if abs(den) > 1e-12*lr*ls                                # not parallel: one point, if it is on both
        t = (w[1]*s[2] - w[2]*s[1])/den
        u = (w[1]*r[2] - w[2]*r[1])/den
        -tol/lr ≤ t ≤ 1 + tol/lr && -tol/ls ≤ u ≤ 1 + tol/ls && push!(ts, clamp(t, 0.0, 1.0))
    elseif abs(w[1]*r[2] - w[2]*r[1]) ≤ tol*lr               # collinear: the other's ends that lie on it
        for q ∈ (c, d)
            t = ((q[1] - a[1])*r[1] + (q[2] - a[2])*r[2])/lr^2
            -tol/lr ≤ t ≤ 1 + tol/lr && push!(ts, clamp(t, 0.0, 1.0))
        end
    end
    return ts
end

"""
    cut(a, bs, scale) -> outlines, or nothing

The part of the zone with outlines `a` (even–odd) outside every zone of `bs` (each a vector of outlines), as
outlines: empty when nothing is left, `nothing` when no zone of `bs` takes any of it. Its vertices are theirs,
exactly, and the points where their outlines cross. `scale` [m], the domain's size, sets the tolerances: points
closer than 1e-9 scale are one, the sides of an edge are probed 1e-7 scale away.
"""
function cut(a, bs, scale)
    tol, ε = 1e-9scale, 1e-7scale
    A, Bs = Polygon(a), Polygon.(bs)                          # their inside rule, `case_points`'s
    segs, own, pts = Tuple{Vector{Float64},Vector{Float64}}[], Bool[], Vector{Float64}[]
    for (k, rs) ∈ enumerate([[a]; bs]), r ∈ rs, i ∈ eachindex(r)
        p, q = r[i], r[mod1(i + 1, length(r))]
        dist(p, q) > tol || continue
        push!(segs, (p, q)); push!(own, k == 1)
        point_index!(pts, p, tol)                             # the vertices first, exactly as they are
    end
    # every edge split where the others meet it: the pieces between consecutive points, each once (both directions
    # are one piece), and whether it lies on the zone's own outline
    pieces = Dict{Tuple{Int,Int},Bool}()
    for (k, (p, q)) ∈ enumerate(segs)
        ts = [0.0, 1.0]
        for (l, (c, d)) ∈ enumerate(segs)
            l == k || append!(ts, meets(p, q, c, d, tol))
        end
        ids = unique([point_index!(pts, t == 0 ? p : t == 1 ? q : p .+ t .* (q .- p), tol) for t ∈ sort(ts)])
        for m ∈ 1:length(ids)-1
            key = minmax(ids[m], ids[m+1])
            pieces[key] = get(pieces, key, false) || own[k]
        end
    end
    # a piece bounds the difference when it lies on exactly one side of it: kept, directed with the difference on its
    # left. Unchanged when the kept pieces are exactly the zone's own
    directed, changed = Tuple{Int,Int}[], false
    for ((i, j), on_a) ∈ pieces
        p, q = pts[i], pts[j]
        d = q .- p
        n, m = [-d[2], d[1]] ./ hypot(d...), (p .+ q) ./ 2
        left, right = (inside(A, x) && !any(B -> inside(B, x), Bs) for x ∈ (m .+ ε .* n, m .- ε .* n))
        changed |= (left != right) != on_a
        left != right && push!(directed, left ? (i, j) : (j, i))
    end
    changed || return nothing
    # chained into outlines: from each piece's end, the unused piece turning most to the left
    outgoing = Dict{Int,Vector{Int}}()
    for (k, (i, _)) ∈ enumerate(directed)
        push!(get!(outgoing, i, Int[]), k)
    end
    used, rings = falses(length(directed)), Vector{Vector{Vector{Float64}}}()
    for k0 ∈ eachindex(directed)
        used[k0] && continue
        ring, k = Int[], k0
        while true
            used[k] = true
            i, j = directed[k]
            push!(ring, i)
            next = [c for c ∈ get(outgoing, j, Int[]) if !used[c]]
            isempty(next) && break
            k = next[argmax([turn(pts[j] .- pts[i], pts[directed[c][2]] .- pts[j]) for c ∈ next])]
        end
        length(ring) ≥ 3 && push!(rings, [copy(pts[i]) for i ∈ ring])
    end
    return rings
end
