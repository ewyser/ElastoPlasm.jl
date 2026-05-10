# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete BSplineBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export BSplineBasis

struct BSplineBasis <: AbstractBasis end

@inline function t1_ϕ∂ϕ(ξ::T2) where {T2}
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
@inline function t2_ϕ∂ϕ(ξ::T2) where {T2}
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
@inline function t3_ϕ∂ϕ(ξ::T2) where {T2}
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
@inline function t4_ϕ∂ϕ(ξ::T2) where {T2}
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
@inline function ϕ∂ϕ(ξ::T2,node_type::T1,Δx::T2) where {T1,T2}
    Δx⁻¹ = T2(1.0)/Δx
    if node_type == T1(1)
        ϕ,∂ϕ = t1_ϕ∂ϕ(ξ*Δx⁻¹)
    elseif node_type == T1(2)
        ϕ,∂ϕ = t2_ϕ∂ϕ(ξ*Δx⁻¹)
    elseif node_type == T1(3)
        ϕ,∂ϕ = t3_ϕ∂ϕ(ξ*Δx⁻¹)
    elseif node_type == T1(4)
        ϕ,∂ϕ = t4_ϕ∂ϕ(ξ*Δx⁻¹)
    else
        error("Invalid type: $(node_type)")
    end   
    return ϕ,∂ϕ*Δx⁻¹ 
end

@inline function basis(mpts::Point{T1,T2,D,B,E,R}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:OneDimension,B<:BSplineBasis,E,R}
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), (T2(0.0),)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[1,ip]-mesh.x[1,no]),mesh.type[1,no],mesh.prprt.h[1])
        # return convolution of basis function
        return T1(no), T2(ϕξ), (T2(∂ϕξ),)
    end
end
@inline function basis(mpts::Point{T1,T2,D,B,E,R}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:TwoDimension,B<:BSplineBasis,E,R} 
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), (T2(0.0), T2(0.0))
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[1,ip]-mesh.x[1,no]),mesh.type[1,no],mesh.prprt.h[1]) 
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[2,ip]-mesh.x[2,no]),mesh.type[2,no],mesh.prprt.h[2])
        # return convolution of basis function
        return T1(no), T2(ϕξ*ϕη), (T2(∂ϕξ*ϕη),T2(ϕξ*∂ϕη),)
    end
end
@inline function basis(mpts::Point{T1,T2,D,B,E,R}, mesh::Mesh{T1,T2,D}, ip::T1, nn::T1) where {T1,T2,D<:ThreeDimension,B<:BSplineBasis,E,R} 
    no = mpts.p2n[nn,ip]
    if iszero(no) 
        return T1(0), T2(0.0), (T2(0.0), T2(0.0), T2(0.0))
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[1,ip]-mesh.x[1,no]),mesh.type[1,no],mesh.prprt.h[1])
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[2,ip]-mesh.x[2,no]),mesh.type[2,no],mesh.prprt.h[2])
        ϕζ,∂ϕζ = ϕ∂ϕ((mpts.x[3,ip]-mesh.x[3,no]),mesh.type[3,no],mesh.prprt.h[3])
        # return convolution of basis function
        return T1(no), T2(ϕξ*ϕη*ϕζ), (T2(∂ϕξ*ϕη*ϕζ),T2(ϕξ*∂ϕη*ϕζ),T2(ϕξ*ϕη*∂ϕζ),)
    end
end