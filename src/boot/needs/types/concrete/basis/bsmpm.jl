# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete BSplineBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export BSplineBasis

struct BSplineBasis{T,D,NN} <: AbstractBasis
    stencils::NTuple{D,UnitRange{T}}
    #
    function BSplineBasis{T,D}() where {T,D}
        stencils = ntuple(_ -> UnitRange{T}(-1:2), D)
        NN       = prod(length.(stencils))
        return new{T,D,NN}(stencils)
    end
end

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

@inline function eval_basis(mpts::Point{T1,T2,1,E,R}, mesh::Mesh{T1,T2,1}, basis::Basis{T1,1,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N  = T2(0.0), SVector{1,T2}(0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),basis.type[1,no],mesh.prprt.h[1])
        # return convolution of basis function
        N, ∂N  = T2(ϕξ), SVector{1,T2}(∂ϕξ)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,2,E,R}, mesh::Mesh{T1,T2,2}, basis::Basis{T1,2,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N  = T2(0.0), SVector{2,T2}(0.0, 0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),basis.type[1,no],mesh.prprt.h[1])
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[ip][2]-mesh.x[no][2]),basis.type[2,no],mesh.prprt.h[2])
        # Construct the convolution of basis function and derivatives
        N, ∂N  = T2(ϕξ*ϕη), SVector{2,T2}(∂ϕξ*ϕη, ϕξ*∂ϕη)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,3,E,R}, mesh::Mesh{T1,T2,3}, basis::Basis{T1,3,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N  = T2(0.0), SVector{3,T2}(0.0, 0.0, 0.0)
    else
        ϕξ,∂ϕξ = ϕ∂ϕ((mpts.x[ip][1]-mesh.x[no][1]),basis.type[1,no],mesh.prprt.h[1])
        ϕη,∂ϕη = ϕ∂ϕ((mpts.x[ip][2]-mesh.x[no][2]),basis.type[2,no],mesh.prprt.h[2])
        ϕζ,∂ϕζ = ϕ∂ϕ((mpts.x[ip][3]-mesh.x[no][3]),basis.type[3,no],mesh.prprt.h[3])
        # Construct the convolution of basis function and derivatives
        N, ∂N  = T2(ϕξ*ϕη*ϕζ), SVector{3,T2}(∂ϕξ*ϕη*ϕζ, ϕξ*∂ϕη*ϕζ, ϕξ*ϕη*∂ϕζ)
    end
    return no, N, ∂N
end


@inline function get_ξηζ(xp::SVector{S,T}, xn::Vector{SVector{S,T}}, h::SVector{S,T},p2n::SVector{NN,Int}) where {S,T,NN}
    nno = ntuple(i -> p2n[i], Val(NN))
    ξηζ = ntuple(Val(NN)) do i
        no = nno[i]
        (
            ntuple(j -> (xp[j] - xn[no][j]) / h[j], Val(S))
        )
    end # ξηζ = ((ξ₁, η₁, ζ₁), (ξ₂, η₂, ζ₂), ..., (ξn, ηn, ζn))
    return nno, ξηζ
end