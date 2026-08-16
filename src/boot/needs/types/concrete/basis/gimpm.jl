# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete LinearBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export GimpBasis

struct GimpBasis <: AbstractBasis end

@inline function S∂S(δx::T2,h::T2,lp::T2) where {T2}                                                       
    S,∂S = T2(0.0),T2(0.0)
    if abs(δx) < lp                       
        S  = T2(1.0)-((T2(4.0)*δx^2+(T2(2.0)*lp)^2)/(T2(8.0)*h*lp))                                   
        ∂S = -((T2(8.0)*δx)/(T2(8.0)*h*lp))                                     
    elseif (abs(δx)>=lp ) && (abs(δx)<=(h-lp))
        S  = T2(1.0)-(abs(δx)/h)                                                       
        ∂S = sign(δx)*(-T2(1.0)/h)                                                   
    elseif (abs(δx)>=(h-lp)) && (abs(δx)< (h+lp))
        S  = ((h+lp-abs(δx))^2)/(T2(4.0)*h*lp)                                       
        ∂S = -sign(δx)*(h+lp-abs(δx))/(T2(2.0)*h*lp)
    end
    return S,∂S    
end
@inline function gimp(mpts::Point{T1,T2,1}, mesh::Mesh{T1,T2,1}, ip::T1, nn::T1) where {T1,T2}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[1,ip]) 
        # convolution of basis function
        N      =  Sx                                     
        ∂Nx    = dSx                                 
        return T1(no), T2(N), T2(∂Nx)
    end
end
@inline function gimp(mpts::Point{T1,T2,2}, mesh::Mesh{T1,T2,2}, ip::T1, nn::T1) where {T1,T2}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0), T2(0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        η      = (mpts.x[ip][2]-mesh.x[no][2])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[ip][1])
        Sy,dSy = S∂S(η,mesh.prprt.h[2],mpts.ℓ[ip][2])
        # convolution of basis function
        N      =  Sx*  Sy                                                                                
        ∂Nx    = dSx*  Sy                                                                                
        ∂Ny    =  Sx* dSy  
        return T1(no), T2(N), T2(∂Nx), T2(∂Ny)
    end
end
@inline function gimp(mpts::Point{T1,T2,3}, mesh::Mesh{T1,T2,3}, ip::T1, nn::T1) where {T1,T2}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0), T2(0.0), T2(0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        η      = (mpts.x[ip][2]-mesh.x[no][2])
        ζ      = (mpts.x[ip][3]-mesh.x[no][3])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[ip][1])
        Sy,dSy = S∂S(η,mesh.prprt.h[2],mpts.ℓ[ip][2])
        Sz,dSz = S∂S(ζ,mesh.prprt.h[3],mpts.ℓ[ip][3])
        # convolution of basis function
        N      =  Sx*  Sy*  Sz                                                                                
        ∂Nx    = dSx*  Sy*  Sz                                                                                
        ∂Ny    =  Sx* dSy*  Sz                                   
        ∂Nz    =  Sx*  Sy* dSz
        return T1(no), T2(N), T2(∂Nx), T2(∂Ny), T2(∂Nz)
    end
end
@inline Base.@propagate_inbounds (::GimpBasis)(mpts, mesh, p, nn) = gimp(mpts, mesh, p, nn)