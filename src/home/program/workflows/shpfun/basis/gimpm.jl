function S∂S(δx,h,lp)                                                         
    S,∂S = 0.0,0.0
    if abs(δx) < lp                       
        S  = 1.0-((4.0*δx^2+(2.0*lp)^2)/(8.0*h*lp))                                   
        ∂S = -((8.0*δx)/(8.0*h*lp))                                     
    elseif (abs(δx)>=   lp ) && (abs(δx)<=(h-lp))
        S  = 1.0-(abs(δx)/h)                                                       
        ∂S = sign(δx)*(-1.0/h)                                                   
    elseif (abs(δx)>=(h-lp)) && (abs(δx)< (h+lp))
        S  = ((h+lp-abs(δx))^2)/(4.0*h*lp)                                       
        ∂S = -sign(δx)*(h+lp-abs(δx))/(2.0*h*lp)
    end
    return S,∂S    
end
@views @kernel inbounds = true function gimpm_1d(mp,mesh)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mp.x[p]-mesh.xn[no]) 
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mp.ℓ[p]) 
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] =  ϕx
            mp.ϕ∂ϕ[nn,p,2] = dϕx
        end
    end
end
@views @kernel inbounds = true function gimpm_2d(mp,mesh)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mp.x[1,p]-mesh.xn[1,no]) 
            η      = (mp.x[2,p]-mesh.xn[2,no])
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mp.ℓ[1,p]) 
            ϕz,dϕz = S∂S(η,mesh.h[2],mp.ℓ[2,p])
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mp.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mp.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz
        end
    end
end
@views @kernel inbounds = true function gimpm_3d(mp,mesh)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mp.nmp
        for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mp.x[1,p]-mesh.xn[1,no]) 
            η      = (mp.x[2,p]-mesh.xn[2,no])
            ζ      = (mp.x[3,p]-mesh.xn[3,no])
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mp.ℓ[1,p])
            ϕy,dϕy = S∂S(η,mesh.h[2],mp.ℓ[2,p])
            ϕz,dϕz = S∂S(ζ,mesh.h[3],mp.ℓ[3,p])
            # convolution of basis function
            mp.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mp.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mp.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mp.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz    
        end
    end
end