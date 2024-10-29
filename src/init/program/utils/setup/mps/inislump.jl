function inislump(meD,cmParam,ni,instr)
    @info "init slump geometry"
    coh0,cohr,phi0,phir,rho0 = cmParam[:c0],cmParam[:cr],cmParam[:ϕ0],cmParam[:ϕr],cmParam[:ρ0]

    lz = 12.80
    wl = 0.15*lz
    if meD.nD == 2
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
        npx,npz     = length(xL),length(zL)
        xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL ))         
        if instr[:GRF]
            c = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])
        else 
            c = ones(size(xp)).*coh0 
        end
        xp,zp,c     = vec(xp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*meD.L[1],a.*x
        xlt,zlt,clt = Float64[],Float64[],Float64[]
        pos         = Float64 
        for mp ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx,Δz = xp[mp]-x[p],zp[mp]-z[p]
                nx,nz = a,-1.0
                if (Δx*nx+Δz*nz)>0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mp]<wl 
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(zlt, zp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(clt, c[mp])
            end
        end
    elseif meD.nD == 3
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        yL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:meD.xB[4]
        zL          = meD.xB[5]+(0.5*meD.h[3]/ni):meD.h[3]/ni:lz-0.5*meD.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        if instr[:GRF]
            c = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])
        else 
            c = ones(size(xp)).*coh0 
        end  
        xp,yp,zp,c  = vec(xp),vec(yp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*meD.L[1],a.*x
        xlt,ylt,zlt = Float64[],Float64[],Float64[]
        clt         = Float64[]
        pos         = Float64 
        for mp ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx = xp[mp]-x[p]
                Δz = zp[mp]-z[p]
                nx = a
                nz = -1.0
                s  = Δx*nx+Δz*nz        
                if s>0.0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mp]<wl 
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mp]) 
                push!(ylt, yp[mp]) 
                push!(zlt, zp[mp]) 
                push!(clt, c[mp])
            end
        end
    end
    
    if meD.nD == 2 
        xp = hcat(xlt,zlt) 
    elseif meD.nD == 3 
        xp = hcat(xlt,ylt,zlt) 
    end
    id = shuffle(collect(1:size(xp,1)))
    coh0   = clt
    cohr   = ones(size(xp,1)).*cohr
    phi    = ones(size(xp,1)).*phi0
    phi[xp[:,end].<=2*wl] .= phir

    return ni,size(xp,1),(;xp=xp,coh0=coh0,cohr=cohr,phi=phi,)
end