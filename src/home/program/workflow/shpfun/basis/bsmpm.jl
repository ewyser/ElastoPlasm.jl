function t1_ϕ∂ϕ(ξ::T2) where {T2}
    ϕ,∂ϕ = T2(0.0),T2(0.0)
    if T2(-2.0)<=ξ<=T2(-1.0) 
        ϕ = T2( 1.0/6.0)     *ξ^3+        ξ^2+T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2( 3.0/6.0)     *ξ^2+T2(2.0)*ξ  +T2(2.0)
    elseif T2(-1.0)<=ξ<=T2(0.0) 
        ϕ = T2(-1.0/6.0)     *ξ^3+        ξ  +T2(1.0)
        ∂ϕ= T2(-3.0/6.0)     *ξ^2+T2(1.0)
    elseif  T2(0.0)<=ξ<=T2(1.0) 
        ϕ = T2( 1.0/6.0)     *ξ^3-        ξ  +T2(1.0)
        ∂ϕ= T2( 3.0/6.0)     *ξ^2-T2(1.0)
    elseif  T2(1.0)<=ξ<=T2(2.0) 
        ϕ = T2(-1.0/6.0)     *ξ^3+        ξ^2-T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2(-3.0/6.0)     *ξ^2+T2(2.0)*ξ  -T2(2.0)
    end   
    return ϕ,∂ϕ
end
function t2_ϕ∂ϕ(ξ::T2) where {T2}
    ϕ,∂ϕ = T2(0.0),T2(0.0)
    if T2(-1.0)<=ξ<=T2(0.0)
        ϕ = T2(-1.0/3.0)     *ξ^3-T2(ξ^2)      +T2(2.0/3.0)
        ∂ϕ= T2(-3.0/3.0)     *ξ^2-T2(2.0)     *ξ
    elseif T2(0.0)<=ξ<=T2(1.0)
        ϕ = T2( 1.0/2.0)     *ξ^3-T2(ξ^2)      +T2(2.0/3.0)
        ∂ϕ= T2( 3.0/2.0)     *ξ^2-T2(2.0)    *ξ
    elseif T2(1.0)<=ξ<=T2(2.0)
        ϕ = T2(-1.0/6.0)     *ξ^3+ξ^2-T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2(-3.0/6.0)     *ξ^2+T2(2.0)    *ξ-T2(2.0)
    end
    return ϕ,∂ϕ
end
function t3_ϕ∂ϕ(ξ::T2) where {T2}
    ϕ,∂ϕ = T2(0.0),T2(0.0)
    if T2(-2.0)<=ξ<=T2(-1.0) 
        ϕ = T2( 1.0/6.0)     *ξ^3+ξ^2      +T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2( 3.0/6.0)     *ξ^2+T2(2.0)*ξ+T2(2.0)
    elseif T2(-1.0)<=ξ<=T2(0.0) 
        ϕ = T2(-1.0/2.0)     *ξ^3-ξ^2      +T2(2.0/3.0)
        ∂ϕ= T2(-3.0/2.0)     *ξ^2-T2(2.0)*ξ
    elseif T2(0.0)<=ξ<=T2(1.0)
        ϕ = T2( 1.0/2.0)     *ξ^3-T2(ξ^2)  +T2(2.0/3.0)
        ∂ϕ= T2( 3.0/2.0)     *ξ^2-T2(2.0)*ξ
    elseif  T2(1.0)<=ξ<=T2(2.0)
        ϕ = T2(-1.0/6.0)     *ξ^3+ξ^2      -T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2(-3.0/6.0)     *ξ^2+T2(2.0)*ξ-T2(2.0)
    end
    return ϕ,∂ϕ
end
function t4_ϕ∂ϕ(ξ::T2) where {T2}
    ϕ,∂ϕ = T2(0.0),T2(0.0)
    if T2(-2.0)<=ξ<=T2(-1.0)
        ϕ = T2( 1.0/6.0)     *ξ^3+ξ^2      +T2(2.0)*ξ+T2(4.0/3.0)
        ∂ϕ= T2( 3.0/6.0)     *ξ^2+T2(2.0)*ξ+T2(2.0)
    elseif T2(-1.0)<=ξ<=T2(0.0)
        ϕ = T2(-1.0/2.0)     *ξ^3-ξ^2      +T2(2.0/3.0)
        ∂ϕ= T2(-3.0/2.0)     *ξ^2-T2(2.0)*ξ
    elseif T2(0.0)<=ξ<=T2(1.0)
        ϕ = T2( 1.0/3.0)     *ξ^3-ξ^2      +T2(2.0/3.0)
        ∂ϕ= T2( 3.0/3.0)     *ξ^2-T2(2.0)*ξ
    end
    return ϕ,∂ϕ
end
function ϕ∂ϕ(ξ::T2,xn::T2,xB::SubArray{T2},Δx::T2) where {T2}
    if xn == xB[1]
        ϕ,∂ϕ = t1_ϕ∂ϕ(ξ)
    elseif xn == (xB[1]+Δx)
        ϕ,∂ϕ = t2_ϕ∂ϕ(ξ)
    elseif (xB[1]+Δx) < xn < (xB[2]-Δx) 
        ϕ,∂ϕ = t3_ϕ∂ϕ(ξ)
    elseif xn == (xB[2]-Δx)
        ϕ,∂ϕ = t4_ϕ∂ϕ(ξ)
    elseif xn==xB[2] 
        ϕ,∂ϕ = t1_ϕ∂ϕ(ξ)
    else
        error("Invalid position: $(xn) not in $(xB)")
    end   
    return ϕ,∂ϕ/Δx 
end
@views @kernel inbounds = true function bsmpm_1d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[p]-mesh.x[no])/mesh.h[1]
            ϕx,dϕx = ϕ∂ϕ(ξ,mesh.x[no],mesh.xB[1:2],mesh.h[1])
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx
            mpts.ϕ∂ϕ[nn,p,2] = dϕx
        end
    end
end
@views @kernel inbounds = true function bsmpm_2d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[1,p]-mesh.x[1,no])/mesh.h[1]
            η      = (mpts.x[2,p]-mesh.x[2,no])/mesh.h[2]
            ϕx,dϕx = ϕ∂ϕ(ξ,mesh.x[1,no],mesh.xB[1,:],mesh.h[1])
            ϕz,dϕz = ϕ∂ϕ(η,mesh.x[2,no],mesh.xB[2,:],mesh.h[2])
            #println("$(typeof(ξ)),$(typeof(η)),$(typeof(ϕx)),$(typeof(ϕz)),$(typeof(dϕx)),$(typeof(dϕz))")
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz   
        end
    end
end
@views @kernel inbounds = true function bsmpm_3d(mpts::Point{T1,T2},mesh::Mesh{T1,T2}) where {T1,T2}
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpts.nmp
        for nn ∈ 1:mesh.nn
            no = mpts.p2n[nn,p]
            if iszero(no) continue end
            # compute basis functions
            ξ      = (mpts.x[1,p]-mesh.x[1,no])/mesh.h[1]
            η      = (mpts.x[2,p]-mesh.x[2,no])/mesh.h[2]
            ζ      = (mpts.x[3,p]-mesh.x[3,no])/mesh.h[3]
            ϕx,dϕx = ϕ∂ϕ(ξ,mesh.x[1,no],mesh.xB[1,:],mesh.h[1])
            ϕy,dϕy = ϕ∂ϕ(η,mesh.x[2,no],mesh.xB[2,:],mesh.h[2])
            ϕz,dϕz = ϕ∂ϕ(ζ,mesh.x[3,no],mesh.xB[3,:],mesh.h[3])
            # convolution of basis function
            mpts.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpts.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpts.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz
        end
    end
end