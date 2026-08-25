"""
    get_slump(mesh::Mesh{T1,T2,D}, mat, solver::S; ni=2, lz=12.80) where {D,S<:AbstractSolver}

Initialize geometry and material point fields for a slump test problem.

# Arguments
- `mesh::Mesh{T1,T2,D}`: Mesh object with geometry and boundary self.
- `mat`: Material parameters (Dict or NamedTuple).
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
    id  = findall(x -> x ≤ lz-(0.5*props.h[end]/ni), out.x[end,:])
    if D == 2
        xp,zp,c     = out.x[1,id],out.x[2,id],out.c0[id]
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*props.L[1],a.*x
        xlt,zlt,clt = Float64[],Float64[],Float64[]
        pos         = Float64
        for mpts ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx,Δz = xp[mpts]-x[p],zp[mpts]-z[p]
                nx,nz = a,-1.0
                if (Δx*nx+Δz*nz)>0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mpts]<wl
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mpts]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(zlt, zp[mpts]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(clt, c[mpts])
            end
        end
    elseif D == 3
        xp,yp,zp,c  = out.x[1,id],out.x[2,id],out.x[3,id],out.c0[id]
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*props.L[1],a.*x
        xlt,ylt,zlt = Float64[],Float64[],Float64[]
        clt         = Float64[]
        pos         = Float64
        for mpts ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx = xp[mpts]-x[p]
                Δz = zp[mpts]-z[p]
                nx = a
                nz = -1.0
                s  = Δx*nx+Δz*nz
                if s>0.0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mpts]<wl
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mpts])
                push!(ylt, yp[mpts])
                push!(zlt, zp[mpts])
                push!(clt, c[mpts])
            end
        end
    end

    if D == 2
        xp = vcat(xlt',zlt')
    elseif D == 3
        xp = vcat(xlt',ylt',zlt')
    end
    nmp    = size(xp,2)
    id     = shuffle(collect(1:nmp))
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