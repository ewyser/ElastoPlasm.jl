#="""
    mpts_populate(mesh, cmpr, instr; ni=2)

Initialize material point fields and coordinates for a mesh, for use in MPM/ElastoPlasm simulations.

# Arguments
- `mesh`: Mesh object containing geometry and boundary information.
- `cmpr`: Material parameters (Dict or NamedTuple), must include at least `:c0` and `:cr`.
- `instr`: Instruction NamedTuple or Dict, may include Gaussian Random Field (GRF) options under `:grf`.
- `ni`: Number of intervals per element (default: 2).

# Returns
- `NamedTuple`: Contains
    - `:ni`   — Number of intervals per element.
    - `:x`    — Matrix of material point coordinates (N×2 for 2D, N×3 for 3D).
    - `:c0`   — Cohesion field (vector, possibly spatially varying if GRF is enabled).
    - (other fields may be added in the future)

# Example
```julia
fields = mpts_populate(mesh, cmpr, instr; ni=4)
@show fields.x
@show fields.c0
```
"""=#
function mpts_populate(mesh,cmpr,instr; ni = 2,)
    out = Dict{Symbol, Any}(:ni => ni)
    if mesh.dim == 2
        x       = collect(mesh.xB[1,1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[1,2])
        z       = collect(mesh.xB[2,1]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[2,2])
        nmp     = [length(x),length(z),length(x)*length(z)]
        xp      = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp      = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])
        out[:x] = hcat(vec(xp),vec(zp))
    elseif mesh.dim == 3
        x       = collect(mesh.xB[1,1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[1,2])
        y       = collect(mesh.xB[2,1]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[2,2])
        z       = collect(mesh.xB[3,1]+(0.5*mesh.h[3]/ni):mesh.h[3]/ni:mesh.xB[3,2])
        nmp     = [length(x),length(y),length(z),length(x)*length(y)*length(z)]
        xp      = repeat(reshape(x,1     ,nmp[1],1     ),nmp[3],1     ,nmp[2])
        yp      = repeat(reshape(y,1     ,1     ,nmp[2]),nmp[3],nmp[1],1     )
        zp      = repeat(reshape(z,nmp[3],1     ,1     ),1     ,nmp[1],nmp[2])
        out[:x] = hcat(vec(xp),vec(yp),vec(zp))
    end
    if instr[:grf][:status]
        if instr[:grf][:covariance] == "gaussian"
            out[:c0] = vec(GRFS_gauss(xp,cmpr[:c0],cmpr[:cr],ni,mesh.h[1]))
        end
        if instr[:grf][:covariance] == "exponential"

        end
    else 
        out[:c0] = ones(size(xp)).*cmpr[:c0]
    end
    return (; (k => v for (k, v) ∈ out)...)
end