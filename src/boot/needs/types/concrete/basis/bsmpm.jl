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

@inline function basis(mpts::Point{T1,T2,1,NN,B,E,R}, mesh::Mesh{T1,T2,1}, ip::T1, nn::T1) where {T1,T2,NN,B<:BSplineBasis,E,R}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        N, ∂N  = T2(0.0), SVector{1,T2}(0.0, 0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),mesh.type[1,no],mesh.prprt.h[1])
        # return convolution of basis function
        N, ∂N  = T2(ϕξ*ϕη), SVector{1,T2}(∂ϕξ*ϕη)
    end
    return no, N, ∂N
end
@inline function basis(mpts::Point{T1,T2,2,NN,B}, mesh::Mesh{T1,T2,2}, ip::T1, nn::T1) where {T1,T2,NN,B<:BSplineBasis}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        N, ∂N  = T2(0.0), SVector{2,T2}(0.0, 0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),mesh.type[1,no],mesh.prprt.h[1])
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[ip][2]-mesh.x[no][2]),mesh.type[2,no],mesh.prprt.h[2])
        # Construct the convolution of basis function and derivatives
        N, ∂N  = T2(ϕξ*ϕη), SVector{2,T2}(∂ϕξ*ϕη, ϕξ*∂ϕη)
    end
    return no, N, ∂N
end
@inline function basis(mpts::Point{T1,T2,3,NN,B,E,R}, mesh::Mesh{T1,T2,3}, ip::T1, nn::T1) where {T1,T2,NN,B<:BSplineBasis,E,R}
    no = mpts.p2n[ip][nn]
    if iszero(no) 
        N, ∂N  = T2(0.0), SVector{3,T2}(0.0, 0.0, 0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),mesh.type[1,no],mesh.prprt.h[1])
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[ip][2]-mesh.x[no][2]),mesh.type[2,no],mesh.prprt.h[2])
        ϕζ,∂ϕζ = ϕ∂ϕ((mpts.x[ip][3]-mesh.x[no][3]),mesh.type[3,no],mesh.prprt.h[3])
        # Construct the convolution of basis function and derivatives
        N, ∂N  = T2(ϕξ*ϕη*ϕζ), SVector{3,T2}(∂ϕξ*ϕη*ϕζ, ϕξ*∂ϕη*ϕζ, ϕξ*ϕη*∂ϕζ)
    end
    return no, N, ∂N
end


@inline function get_ξηζ(xp::SVector{S,T}, xn::Vector{SVector{S,T}}, h::SVector{S,T},p2n::SVector{N,Int}) where {S,T,N}
    nno = ntuple(i -> p2n[i], Val(N))
    ξηζ = ntuple(Val(N)) do i
        no = nno[i]
        (
            ntuple(j -> (xp[j] - xn[no][j]) / h[j], Val(S))
        )
    end # ξηζ = ((ξ₁, η₁, ζ₁), (ξ₂, η₂, ζ₂), ..., (ξn, ηn, ζn))
    return nno, ξηζ
end