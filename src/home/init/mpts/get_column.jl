"""
    get_column(mesh::Mesh{T1,T2,D}, mat, ni; ℓ₀=0.0) where {D}

Initialize geometry and material point fields for an elastic self-weight column test.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object with geometry and boundary self.
- `mat`: Material parameters (Dict or NamedTuple).
- `ni`: Number of intervals per element.
- `ℓ₀`: Optional, column height (default: 0.0).

# Returns
- `ni`: Number of intervals per element.
- `nmp`: Number of material points.
- `fields`: NamedTuple with coordinates and material properties.
"""
function get_column(mesh::Mesh{T1,T2,D},mat,ni;ℓ₀=0.0) where {T1,T2,D}
    @info "Init elastic collumn geometry"
    if D == 2
        x          = collect(mesh.prprt.xB[1,1]+(0.5*mesh.prprt.h[1]/ni):mesh.prprt.h[1]/ni:mesh.prprt.xB[1,2])
        z          = collect(mesh.prprt.xB[2,1]+(0.5*mesh.prprt.h[2]/ni):mesh.prprt.h[2]/ni:ℓ₀        )
        nmp        = [length(x),length(z),length(x)*length(z)]
        xp         = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp         = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])
    elseif D == 3
        # TODO: 3D branch unimplemented — xp/yp/zp never assigned, would throw UndefVarError
        #=
        xL          = mesh.prprt.xB[1,1]+(0.5*mesh.prprt.h[1]/ni):mesh.prprt.h[1]/ni:mesh.prprt.xB[1,2]
        yL          = mesh.prprt.xB[2,1]+(0.5*mesh.prprt.h[2]/ni):mesh.prprt.h[2]/ni:mesh.prprt.xB[2,2]
        zL          = mesh.prprt.xB[3,1]+(0.5*mesh.prprt.h[3]/ni):mesh.prprt.h[3]/ni:ℓ₀-0.5*mesh.prprt.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        xp,yp,zp    = vec(xp),vec(yp),vec(zp)=#
    end
    if D == 2
        xp = vcat(vec(xp)',vec(zp)')
    elseif D == 3
        xp = vcat(vec(xp)',vec(yp)',vec(zp)')
    end
    nmp  = size(xp,2)
    return (; xp, ni, nmp, material_fields(mat, nmp)...)
end