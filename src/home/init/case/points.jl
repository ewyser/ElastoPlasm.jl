export case_points

"""
    case_points(case, mesh, solver) -> PointData

The material points of a case: the candidate points of `mpts_populate` (`case.ni`×`case.ni` per
element, at sub-cell centres) that lie in a zone, each taking the last zone containing it. Every
constant is the case's material (`case.material`) unless the point's zone sets it; the cohesion
of points whose zone does not set `c0` is the random field of `solver.grf` when enabled
(`property_field`, bounded below by the material's `cr`), evaluated over all points at once. A
zone's velocity `v` has `D` components, or two (x, z) in 3D.
"""
function case_points(case::Case, mesh::Mesh{T1,T2,D}, solver) where {T1,T2,D}
    cand = mpts_populate(mesh.prprt; ni=case.ni)
    zone = [something(findlast(z -> inside(z.shape, view(cand, :, p)), case.zones), 0) for p ∈ axes(cand, 2)]
    keep = zone .> 0
    xp, zone = cand[:, keep], zone[keep]
    nmp, m = size(xp, 2), case.material
    f = Dict{Symbol,Vector{Float64}}(k => fill(getfield(m, k), nmp) for k ∈ keys(m))
    f[:c0] = property_field(xp, m.c0, solver.grf; floor=m.cr)
    vp = zeros(D, nmp)
    for (i, z) ∈ enumerate(case.zones), (k, v) ∈ z.props
        id = zone .== i
        if k == :v                                          # (x, z) in 3D: no velocity along y, as 2D shapes extrude
            w = length(v) == D ? v : D == 3 && length(v) == 2 ? [v[1], 0.0, v[2]] :
                throw(ArgumentError("zone \"$(z.name)\": v needs $D (or, in 3D, 2: x, z) components"))
            vp[:, id] .= w
        else
            f[k][id] .= v
        end
    end
    return PointData(xp, vp, case.ni, zone, f[:c0], f[:cr], f[:phi], f[:Hp], f[:E], f[:nu], f[:rho0], f[:n0], f[:T], f[:cp], f[:k])
end
