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
        return ϕ,∂ϕ/Δx  
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
        return ϕ,∂ϕ/Δx
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
        return ϕ,∂ϕ/Δx
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
        return ϕ,∂ϕ/Δx
    end    
end
@views @kernel inbounds = true function bsmpm_1d(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1]) 
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.xB[1:2],meD.h[1])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx
            mpD.ϕ∂ϕ[nn,p,2] = dϕx
        end
    end
end
@views @kernel inbounds = true function bsmpm_2d(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1]) 
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.xB[1:2],meD.h[1])
            ϕz,dϕz = ϕ∂ϕ(η/meD.h[2],meD.xn[no,2],meD.xB[3:4],meD.h[2])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz   
        end
    end
end
@views @kernel inbounds = true function bsmpm_3d(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1])
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ζ      = (mpD.x[p,3]-meD.xn[no,3])
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.xB[1:2],meD.h[1])
            ϕy,dϕy = ϕ∂ϕ(η/meD.h[2],meD.xn[no,2],meD.xB[3:4],meD.h[2])
            ϕz,dϕz = ϕ∂ϕ(ζ/meD.h[3],meD.xn[no,3],meD.xB[5:6],meD.h[3])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpD.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz
        end
    end
end