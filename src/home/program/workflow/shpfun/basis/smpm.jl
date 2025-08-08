function N∂N(δx,h)                                                       
    if -h<δx<=0.0                       
        N,∂N = 1.0+δx/h,1.0/h                                                                         
    elseif 0.0<δx<=h
        N,∂N = 1.0-δx/h,-1.0/h                                                                                                          
    else
        N,∂N = 0.0,0.0                                                                                                  
    end
    return N,∂N    
end
@views @kernel inbounds = true function smpm_1d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[p]-mesh.x[no])
            ϕx,dϕx = N∂N(ξ,mesh.h[1]    )
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx
            mpts.ϕ∂ϕ[nn,p,2] = dϕx
        end  
    end
end
@views @kernel inbounds = true function smpm_2d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[1,p]-mesh.x[1,no]) 
            η      = (mpts.x[2,p]-mesh.x[2,no])
            ϕx,dϕx = N∂N(ξ,mesh.h[1]        )
            ϕz,dϕz = N∂N(η,mesh.h[2]        )
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz
        end
    end
end
@views @kernel inbounds = true function smpm_3d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
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
            ϕx,dϕx = N∂N(ξ,mesh.h[1]        )
            ϕy,dϕy = N∂N(η,mesh.h[2]        )
            ϕz,dϕz = N∂N(ζ,mesh.h[3]        )
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpts.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz  
        end
    end
end