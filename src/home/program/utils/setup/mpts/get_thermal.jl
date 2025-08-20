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
function get_thermal(mesh,cmp,ni;ℓ₀=0.0)
    @info "Init thermal diffusion"
    coh0,cohr,phi0= cmp[:c0],cmp[:cr],cmp[:ϕ0]
    if mesh.dim == 2
        x          = collect(mesh.xB[1,1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[1,2])
        z          = collect(mesh.xB[2,1]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:ℓ₀        )
        nmp        = [length(x),length(z),length(x)*length(z)]
        xp         = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp         = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])
    elseif mesh.dim == 3
        #=
        xL          = mesh.xB[1,1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[1,2]
        yL          = mesh.xB[2,1]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[2,2]
        zL          = mesh.xB[3,1]+(0.5*mesh.h[3]/ni):mesh.h[3]/ni:ℓ₀-0.5*mesh.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        xp,yp,zp    = vec(xp),vec(yp),vec(zp)=#
    end
    if mesh.dim == 2 
        xp = vcat(vec(xp)',vec(zp)') 
    elseif mesh.dim == 3 
        xp = vcat(vec(xp)',vec(yp)',vec(zp)') 
    end
    nmp  = size(xp,2)
    id   = shuffle(collect(1:nmp))
    coh0 = ones(nmp).*coh0
    cohr = ones(nmp).*cohr
    phi  = ones(nmp).*phi0
    return (;xp=xp,coh0=coh0,cohr=cohr,phi=phi,ni=ni,nmp=nmp)
end