function which(xn,xB,Δx)
    if xn == xB[1] ||  xn==xB[2] 
        return 1::Int64
    elseif xn == (xB[1]+Δx)
        return 2::Int64
    elseif (xB[1]+Δx) < xn < (xB[2]-Δx) 
        return 3::Int64
    elseif xn == (xB[2]-Δx)
        return 4::Int64
    else
        return 0::Int64
    end
end
function ϕ∂ϕ(ξ,xn,xB,Δx)
    ϕ,∂ϕ = 0.0,0.0
    if which(xn,xB,Δx) == 1
        if -2.0<=ξ<=-1.0 
            ϕ = 1.0/6.0     *ξ^3+     ξ^2   +2.0*ξ    +4.0/3.0
            ∂ϕ= 3.0/6.0     *ξ^2+2.0 *ξ     +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/6.0     *ξ^3           +    ξ    +1.0
            ∂ϕ= -3.0/6.0     *ξ^2           +  1.0
        elseif  0.0<=ξ<= 1.0 
            ϕ =  1.0/6.0     *ξ^3           -    ξ    +1.0
            ∂ϕ=  3.0/6.0     *ξ^2           -  1.0
        elseif  1.0<=ξ<= 2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2  -2.0*ξ    +4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ    -2.0
        end  
    elseif which(xn,xB,Δx) == 2
        if -1.0<=ξ<=0.0 
            ϕ = -1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/3.0     *ξ^2-2.0 *ξ
        elseif 0.0<=ξ<=1.0 
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif 1.0<=ξ<=2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif which(xn,xB,Δx) == 3
        if -2.0<=ξ<=-1.0 
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ
        elseif  0.0<=ξ<=1.0
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif  1.0<=ξ<=2.0    
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif which(xn,xB,Δx) == 4
        if -2.0<=ξ<=-1.0
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0 
        elseif -1.0<=ξ<=0.0
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ      
        elseif 0.0<=ξ<=1.0
            ϕ =  1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/3.0     *ξ^2-2.0 *ξ      
        end
    end   
    return ϕ,∂ϕ/Δx 
end
@views @kernel inbounds = true function bsmpm_1d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute basis functions
            ξ      = (mp.x[p]-mesh.x[no]) 
            ϕx,dϕx = ϕ∂ϕ(ξ/mesh.h[1],mesh.x[no],mesh.xB[1:2],mesh.h[1])
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] =  ϕx
            mp.ϕ∂ϕ[nn,p,2] = dϕx
        end
    end
end
@views @kernel inbounds = true function bsmpm_2d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute basis functions
            ξ      = T2(mp.x[1,p]-mesh.x[1,no]) 
            η      = T2(mp.x[2,p]-mesh.x[2,no])
            ϕx,dϕx = ϕ∂ϕ(ξ/mesh.h[1],mesh.x[1,no],mesh.xB[1:2],mesh.h[1])
            ϕz,dϕz = ϕ∂ϕ(η/mesh.h[2],mesh.x[2,no],mesh.xB[3:4],mesh.h[2])
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] = T2( ϕx*  ϕz)                                        
            mp.ϕ∂ϕ[nn,p,2] = T2(dϕx*  ϕz)                                        
            mp.ϕ∂ϕ[nn,p,3] = T2( ϕx* dϕz)   
        end
    end
end
@views @kernel inbounds = true function bsmpm_3d(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for nn ∈ 1:mesh.nn
            no = mp.p2n[nn,p]
            if no < 1 continue end
            # compute basis functions
            ξ      = (mp.x[1,p]-mesh.x[1,no]) 
            η      = (mp.x[2,p]-mesh.x[2,no])
            ζ      = (mp.x[3,p]-mesh.x[3,no])
            ϕx,dϕx = ϕ∂ϕ(ξ/mesh.h[1],mesh.x[1,no],mesh.xB[1:2],mesh.h[1])
            ϕy,dϕy = ϕ∂ϕ(η/mesh.h[2],mesh.x[2,no],mesh.xB[3:4],mesh.h[2])
            ϕz,dϕz = ϕ∂ϕ(ζ/mesh.h[3],mesh.x[3,no],mesh.xB[5:6],mesh.h[3])
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mp.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mp.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mp.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz
        end
    end
end