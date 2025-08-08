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
@views @kernel inbounds = true function gimpm_1d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[p]-mesh.x[no]) 
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mpts.ℓ[p]) 
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx
            mpts.ϕ∂ϕ[nn,p,2] = dϕx
        end
    end
end
@views @kernel inbounds = true function gimpm_2d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[1,p]-mesh.x[1,no]) 
            η      = (mpts.x[2,p]-mesh.x[2,no])
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mpts.ℓ[1,p]) 
            ϕz,dϕz = S∂S(η,mesh.h[2],mpts.ℓ[2,p])
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz
        end
    end
end
@views @kernel inbounds = true function gimpm_3d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[1,p]-mesh.x[1,no]) 
            η      = (mpts.x[2,p]-mesh.x[2,no])
            ζ      = (mpts.x[3,p]-mesh.x[3,no])
            ϕx,dϕx = S∂S(ξ,mesh.h[1],mpts.ℓ[1,p])
            ϕy,dϕy = S∂S(η,mesh.h[2],mpts.ℓ[2,p])
            ϕz,dϕz = S∂S(ζ,mesh.h[3],mpts.ℓ[3,p])
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpts.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz    
        end
    end
end