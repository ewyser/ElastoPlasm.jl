# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete LinearBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export LinearBasis

struct LinearBasis <: AbstractBasis end

@inline function N∂N(δx::T2,h::T2) where {T2}                                                     
    if -h < δx <= T2(0.0)                       
        N,∂N = T2(1.0)+δx/h,T2(1.0)/h                                                                         
    elseif T2(0.0) < δx <= h
        N,∂N = T2(1.0)-δx/h,T2(-1.0)/h                                                                                                          
    else
        N,∂N = T2(0.0),T2(0.0)                                                                                                  
    end
    return N,∂N    
end


@inline function linear(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:OneDimension}
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0)
    else
        φξ,∂φξ = N∂N((mpts.x[ip][1]-mesh.x[no][1]),mesh.prprt.h[1])
        # return convolution of basis function
        return T1(no), T2(ϕξ), (T2(∂ϕξ),)
    end
end
@inline function linear(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:TwoDimension}
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0), T2(0.0)
    else
        φξ,∂φξ = N∂N((mpts.x[1,ip]-mesh.x[no][1]),mesh.prprt.h[1])
        φη,∂φη = N∂N((mpts.x[2,ip]-mesh.x[no][2]),mesh.prprt.h[2])
        # return convolution of basis function
        return T1(no), T2(ϕξ*ϕη), (T2(∂ϕξ*ϕη),T2(ϕξ*∂ϕη),)
    end
end
@inline function linear(mpts::Point{T1,T2,D}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:ThreeDimension}
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), T2(0.0), T2(0.0), T2(0.0)
    else
        φξ,∂φξ = N∂N((mpts.x[1,ip]-mesh.x[no][1]),mesh.prprt.h[1])
        φη,∂φη = N∂N((mpts.x[2,ip]-mesh.x[no][2]),mesh.prprt.h[2])
        φζ,∂φζ = N∂N((mpts.x[3,ip]-mesh.x[no][3]),mesh.prprt.h[3])
        # return convolution of basis function
        return T1(no), T2(ϕξ*ϕη*ϕζ), (T2(∂ϕξ*ϕη*ϕζ),T2(ϕξ*∂ϕη*ϕζ),T2(ϕξ*ϕη*∂ϕζ),)
    end
end
@inline Base.@propagate_inbounds (::LinearBasis)(mpts, mesh, p, nn) = linear(mpts, mesh, p, nn)