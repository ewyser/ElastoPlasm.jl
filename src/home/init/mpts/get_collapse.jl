"""
    get_collapse(mesh::Mesh{T1,T2,D}, mat, solver::S; ni=2, w0=0.0, h0=0.0, ϕ=mat[:ϕ0], c0=mat[:c0]) where {D}

Initialize geometry and material point fields for a dry granular block-collapse
problem: a plain rectangular block of width `w0` and height `h0`, placed at the
domain's lower-left corner, released under gravity with no pre-existing slope —
mirrors MaterialPointSolver.jl's `2d_collapse.jl` reference scenario. Contrast with
`get_slump` (pre-carved sloped geometry) and `get_column` (1-D elastic self-weight
column, no plasticity).

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object with geometry and boundary info.
- `mat`: Material parameters (Dict or NamedTuple), see `setup_material_constants`.
- `solver::S`: Solver instance; may enable GRF cohesion via `solver.grf` (see `mpts_populate`).
- `ni`: Number of intervals per element (default: 2).
- `w0`,`h0`: Initial block width/height — required, no generic default (block size is
  scenario-specific); must be ≤ the mesh domain size.
- `ϕ`: Friction angle [rad] (default: `mat[:ϕ0]`) — overridable since
  `setup_material_constants` doesn't expose cohesion/friction as kwargs.
- `c0`: Cohesion [Pa] (default: `mat[:c0]`) — overridable, same reason; `c0=0.0` gives
  a cohesionless dry-sand block, matching MaterialPointSolver.jl's own comparison case.

# Returns
- `NamedTuple`: `(;xp, coh0, cohr, phi, T, c, k, ni, nmp)`, matching `get_slump`/
  `get_column`'s shape (for `setup_mpts` compatibility).
"""
function get_collapse(mesh::Mesh{T1,T2,D},mat,solver::S; ni=2, w0=0.0, h0=0.0, ϕ=mat[:ϕ0], c0=mat[:c0]) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    @info "Init granular block collapse geometry"
    D == 3 && error("get_collapse: 3D not yet implemented (2D only, matches the MaterialPointSolver.jl comparison scenario)")
    out  = mpts_populate(mesh.prprt, mat, solver; ni=ni)
    xrow = out.x[1,:]
    zrow = out.x[2,:]
    id   = findall(i -> (xrow[i] ≤ w0) && (zrow[i] ≤ h0), eachindex(xrow))
    xp   = out.x[:,id]
    nmp  = length(id)

    coh0 = fill(c0,nmp)
    cohr = fill(c0,nmp)
    phi  = fill(ϕ,nmp)

    c    = fill(mat[:specific_heat_capacity],nmp)
    k    = fill(mat[:thermal_conductivity],nmp)
    T    = fill(mat[:initial_temperature],nmp)

    return (;xp=xp,coh0=coh0,cohr=cohr,phi=phi,T=T,c=c,k=k,ni=ni,nmp=nmp)
end
