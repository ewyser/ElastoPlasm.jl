# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete LinearBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export GimpBasis

struct GimpBasis <: AbstractBasis end

stencil_range(::GimpBasis, ::Type{T1}) where {T1} = T1(-1):T1(2)

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
@inline function eval_basis(mpts::Point{T1,T2,1,E,R}, mesh::Mesh{T1,T2,1}, basis::Basis{T1,1,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{1,T2}(0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[ip][1])
        # convolution of basis function
        N, ∂N  = T2(Sx), SVector{1,T2}(dSx)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,2,E,R}, mesh::Mesh{T1,T2,2}, basis::Basis{T1,2,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{2,T2}(0.0, 0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        η      = (mpts.x[ip][2]-mesh.x[no][2])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[ip][1])
        Sy,dSy = S∂S(η,mesh.prprt.h[2],mpts.ℓ[ip][2])
        # convolution of basis function
        N, ∂N  = T2(Sx*Sy), SVector{2,T2}(dSx*Sy, Sx*dSy)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,3,E,R}, mesh::Mesh{T1,T2,3}, basis::Basis{T1,3,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{3,T2}(0.0, 0.0, 0.0)
    else
        ξ      = (mpts.x[ip][1]-mesh.x[no][1])
        η      = (mpts.x[ip][2]-mesh.x[no][2])
        ζ      = (mpts.x[ip][3]-mesh.x[no][3])
        Sx,dSx = S∂S(ξ,mesh.prprt.h[1],mpts.ℓ[ip][1])
        Sy,dSy = S∂S(η,mesh.prprt.h[2],mpts.ℓ[ip][2])
        Sz,dSz = S∂S(ζ,mesh.prprt.h[3],mpts.ℓ[ip][3])
        # convolution of basis function
        N, ∂N  = T2(Sx*Sy*Sz), SVector{3,T2}(dSx*Sy*Sz, Sx*dSy*Sz, Sx*Sy*dSz)
    end
    return no, N, ∂N
end