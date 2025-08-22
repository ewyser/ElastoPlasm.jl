#="""
    geom_collapse(mesh, cmp, ni; ℓ₀=0.0)

Initialize geometry and material point fields for a column collapse problem.

# Arguments
- `mesh`: Mesh object with geometry and boundary info.
- `cmp`: Material parameters (Dict or NamedTuple).
- `ni`: Number of intervals per element.
- `ℓ₀`: Optional, domain height (default: 0.0).

# Returns
- `ni`: Number of intervals per element.
- `nmp`: Number of material points.
- `fields`: NamedTuple with coordinates and material properties.
"""=#
function get_thermal(mesh,cmpr,instr; ni = 2, )
    props = mesh.prprt
    #@info "Init slump geometry"
    out = mpts_populate(props,cmpr,instr; ni=ni)

    logic1(i) = out.x[end, i] ≤ mesh.prprt.xB[end,2]
    logic2(i) = 0.25*mesh.prprt.xB[1,2] ≤ out.x[1, i] ≤ 0.75*mesh.prprt.xB[1,2]
    id = findall(i -> logic1(i) && logic2(i), axes(out.x,2))
    xp = copy(out.x[:,id])
    #xp = copy(out.x)


    nmp  = size(xp,2)

    coh0 = ones(nmp).*cmpr[:c0]
    cohr = ones(nmp).*cmpr[:cr]
    phi  = ones(nmp).*cmpr[:ϕ0]

    c    = ones(nmp).*cmpr[:specific_heat_capacity]
    k    = ones(nmp).*cmpr[:thermal_conductivity]
    T    = zeros(nmp).*cmpr[:initial_temperature]
    T[xp[end,:].<=0.5*mesh.prprt.xB[end,2]] .= cmpr[:initial_temperature]+60.0

#=
    id = out.x[1, :] .<= mesh.prprt.xB[1,1] + mesh.prprt.h[1]
    T[id] .= 100.0
    id = out.x[1, :] .>= mesh.prprt.xB[1,2] - mesh.prprt.h[1]
    T[id] .= 100.0

    id = out.x[2, :] .<= mesh.prprt.xB[2,1] + mesh.prprt.h[2]
    #T[id] .= 100.0
    id = out.x[2, :] .>= mesh.prprt.xB[2,2] - mesh.prprt.h[2]
    T[id] .= 100.0
=#


    return (;xp=xp,coh0=coh0,cohr=cohr,phi=phi,T=T,c=c,k=k,ni=ni,nmp=nmp)
end