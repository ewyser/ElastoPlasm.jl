function inislump(mesh,cmp,instr; ni = 2, lz = 12.80)
    #@info "Init slump geometry"
    coh0,cohr,phi0,phir,rho0 = cmp[:c0],cmp[:cr],cmp[:ϕ0],cmp[:ϕr],cmp[:ρ0]

    wl = 0.15*lz
    if mesh.dim == 2
        x          = collect(mesh.xB[1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[2])
        z          = collect(mesh.xB[3]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:lz-0.5*mesh.h[2]/ni)
        nmp        = [length(x),length(z),length(x)*length(z)]
        xp         = repeat(reshape(x,1     ,nmp[1]),nmp[2],1     )
        zp         = repeat(reshape(z,nmp[2],1     ),1     ,nmp[1])

        if instr[:grf][:status]
            if instr[:grf][:covariance] == "gaussian"
                c = GRFS_gauss(xp,coh0,cohr,ni,mesh.h[1])
            end
            if instr[:grf][:covariance] == "exponential"

            end
        else 
            c = ones(size(xp)).*coh0 
        end
        xp,zp,c     = vec(xp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*mesh.L[1],a.*x
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
    elseif mesh.dim == 3
        xL          = mesh.xB[1]+(0.5*mesh.h[1]/ni):mesh.h[1]/ni:mesh.xB[2]
        yL          = mesh.xB[3]+(0.5*mesh.h[2]/ni):mesh.h[2]/ni:mesh.xB[4]
        zL          = mesh.xB[5]+(0.5*mesh.h[3]/ni):mesh.h[3]/ni:lz-0.5*mesh.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        if instr[:grf][:status]
            c = GRFS_gauss(xp,coh0,cohr,ni,mesh.h[1])
        else 
            c = ones(size(xp)).*coh0 
        end  
        xp,yp,zp,c  = vec(xp),vec(yp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*mesh.L[1],a.*x
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
    
    if mesh.dim == 2 
        xp = vcat(xlt',zlt') 
    elseif mesh.dim == 3 
        xp = vcat(xlt',ylt',zlt') 
    end
    nmp    = size(xp,2)
    id     = shuffle(collect(1:nmp))
    coh0   = clt
    cohr   = ones(nmp).*cohr
    phi    = ones(nmp).*phi0
    phi[xp[end,:].<=2*wl] .= phir

    return ni,nmp,(;xp=xp,coh0=coh0,cohr=cohr,phi=phi,)
end