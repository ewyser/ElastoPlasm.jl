"""
    get_slump(mesh::Mesh{T1,T2,D}, mat, solver::S; ni=2, lz=12.80) where {D,S<:AbstractSolver}

Initialize geometry and material point fields for a slump test problem.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object with geometry and boundary self.
- `mat`: Material parameters (NamedTuple, see `setup_material_constants`).
- `solver::S`: Solver instance (e.g. `ExplicitSolver`, may include GRF options).
- `ni`: Number of intervals per element (default: 2).
- `lz`: Domain height (default: 12.80).

# Returns
- `ni`: Number of intervals per element.
- `nmp`: Number of material points.
- `fields`: NamedTuple with coordinates and material properties.
"""
function get_slump(mesh::Mesh{T1,T2,D}, mat, solver::S; ni = 2, lz = 12.80) where {T1,T2,D,S<:AbstractSolver}
    props = mesh.prprt
    out = mpts_populate(props,mat,solver; ni=ni)
    wl  = 0.15*lz
    # keep points below the slump height lz (last coordinate is vertical)
    id  = findall(x -> x ≤ lz-(0.5*props.h[end]/ni), out.x[end,:])
    xp,c = out.x[:,id],out.c0[id]
    # slope line z = a·(x - L/2) through (xs, zs); keep points on its inner side, plus the base layer z < wl
    a       = -1.25
    xs,zs   = maximum(xp[1,:])+0.5*props.L[1], a*maximum(xp[1,:])
    keep    = [(xp[1,p]-xs)*a+(xp[end,p]-zs)*(-1.0) > 0 || xp[end,p] < wl for p ∈ axes(xp,2)]
    xp,clt  = xp[:,keep],c[keep]
    nmp    = size(xp,2)
    coh0   = clt
    cohr   = ones(nmp).*mat[:cr]
    phi    = ones(nmp).*mat[:ϕ0]
    phi[xp[end,:].<=2*wl] .= mat[:ϕr]

    c      = ones(nmp).*mat[:specific_heat_capacity]
    k      = ones(nmp).*mat[:thermal_conductivity]
    T      = ones(nmp).*mat[:initial_temperature]
    T[xp[end,:].<=2*wl] .= 3.0*mat[:initial_temperature]

    return (;xp=xp,coh0=coh0,cohr=cohr,phi=phi,T=T,c=c,k=k,ni=ni,nmp=nmp)
end