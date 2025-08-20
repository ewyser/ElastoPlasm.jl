#="""
    get_slump(props, cmp, instr; ni=2, lz=12.80)

Initialize geometry and material point fields for a slump test problem.

# Arguments
- `props`: props object with geometry and boundary info.
- `cmp`: Material parameters (Dict or NamedTuple).
- `instr`: Instruction dictionary (may include GRF options).
- `ni`: Number of intervals per element (default: 2).
- `lz`: Domain height (default: 12.80).

# Returns
- `ni`: Number of intervals per element.
- `nmp`: Number of material points.
- `fields`: NamedTuple with coordinates and material properties.
"""=#
function get_slump(mesh,cmpr,instr; ni = 2, lz = 12.80 )
    props = mesh.prprt
    #@info "Init slump geometry"
    out = mpts_populate(props,cmpr,instr; ni=ni)
    wl  = 0.15*lz
    id  = findall(x -> x ≤ lz-(0.5*props.h[end]/ni), out.x[end,:])
    if props.dim == 2
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
    elseif props.dim == 3
        xp,yp,zp,c  = out.x[id,1],out.x[id,2],out.x[id,3],out.c0[id]
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
    
    if props.dim == 2 
        xp = vcat(xlt',zlt') 
    elseif props.dim == 3 
        xp = vcat(xlt',ylt',zlt') 
    end
    nmp    = size(xp,2)
    id     = shuffle(collect(1:nmp))
    coh0   = clt
    cohr   = ones(nmp).*cmpr[:cr]
    phi    = ones(nmp).*cmpr[:ϕ0]
    phi[xp[end,:].<=2*wl] .= cmpr[:ϕr]


    c      = ones(nmp).*cmpr[:specific_heat_capacity]
    k      = ones(nmp).*cmpr[:thermal_conductivity]
    T      = ones(nmp).*cmpr[:initial_temperature]
    T[xp[end,:].<=2*wl] .= 3.0*cmpr[:initial_temperature]

    return (;xp=xp,coh0=coh0,cohr=cohr,phi=phi,T=T,c=c,k=k,ni=ni,nmp=nmp)
end