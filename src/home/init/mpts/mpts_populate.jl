"""
    mpts_populate(props, cmpr, solver::S; ni=2) where {S<:AbstractSolver}

Initialize candidate material point fields and coordinates on a regular per-element grid, for use in MPM/ElastoPlasm simulations.

# Arguments
- `props`: `MeshProperties` object containing geometry and boundary information.
- `cmpr`: Material parameters (Dict or NamedTuple), must include at least `:c0` and `:cr`.
- `solver::S`: Solver instance (e.g. `ExplicitSolver`), may enable Gaussian Random Field (GRF) cohesion via `solver.grf`.
- `ni`: Number of intervals per element (default: 2).

# Returns
- `NamedTuple`: Contains
    - `:ni`   — Number of intervals per element.
    - `:x`    — Matrix of material point coordinates (N×2 for 2D, N×3 for 3D).
    - `:c0`   — Cohesion field (vector, possibly spatially varying if GRF is enabled).
    - (other fields may be added in the future)

# Example
```julia
fields = mpts_populate(props, cmpr, solver; ni=4)
@show fields.x
@show fields.c0
```
"""
function mpts_populate(props,cmpr,solver::S; ni = 2,) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    out  = Dict{Symbol, Any}(:ni => ni)
    if D == 2
        x       = collect(props.xB[1,1]+(0.5*props.h[1]/ni):props.h[1]/ni:props.xB[1,2]-(0.5*props.h[1]/ni))
        z       = collect(props.xB[2,1]+(0.5*props.h[2]/ni):props.h[2]/ni:props.xB[2,2]-(0.5*props.h[2]/ni))
        nmp     = [length(x),length(z),length(x)*length(z)]
        xp      = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp      = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])
        out[:x] = vcat(vec(xp)',vec(zp)')
    elseif D == 3
        x       = collect(props.xB[1,1]+(0.5*props.h[1]/ni):props.h[1]/ni:props.xB[1,2]-(0.5*props.h[1]/ni))
        y       = collect(props.xB[2,1]+(0.5*props.h[2]/ni):props.h[2]/ni:props.xB[2,2]-(0.5*props.h[2]/ni))
        z       = collect(props.xB[3,1]+(0.5*props.h[3]/ni):props.h[3]/ni:props.xB[3,2]-(0.5*props.h[3]/ni))
        nmp     = [length(x),length(y),length(z),length(x)*length(y)*length(z)]
        xp      = repeat(reshape(x,1     ,nmp[1],1     ),nmp[3],1     ,nmp[2])
        yp      = repeat(reshape(y,1     ,1     ,nmp[2]),nmp[3],nmp[1],1     )
        zp      = repeat(reshape(z,nmp[3],1     ,1     ),1     ,nmp[1],nmp[2])
        out[:x] = vcat(vec(xp)',vec(yp)',vec(zp)')
    end
    if solver.grf.status
        if solver.grf.covariance == "gaussian"
            out[:c0] = vec(GRFS_gauss(xp,cmpr[:c0],cmpr[:cr],ni,props.h[1]))
        end
        if solver.grf.covariance == "exponential"

        end
    else 
        out[:c0] = ones(size(xp)).*cmpr[:c0]
    end
    return (; (k => v for (k, v) ∈ out)...)
end