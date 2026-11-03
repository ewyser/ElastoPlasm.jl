# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete BSplineBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export BSplineBasis

"""
    BSplineBasis{T1,D,NN} <: AbstractBasisKind

Cubic B-spline shape functions, modified near the grid's boundary by the nodes' boundary layer (`basis.type`): a
basis kind, what `Basis`'s `kind` holds.

# Fields
- `stencils::NTuple{D,UnitRange{T1}}`: per axis, the offsets of the nodes a point reaches from its element's first node
  (`-1:2`, 4 nodes per axis).

`T1` is the index type, `D` the dimension, `NN` the number of nodes a point reaches (`prod(length.(stencils))`).

# Constructors
    BSplineBasis{T1,D}()

`get_basis` builds it from the option `basis.which` (`"bsmpm"`).
"""
struct BSplineBasis{T1,D,NN} <: AbstractBasisKind
    stencils::NTuple{D,UnitRange{T1}}
    #
    function BSplineBasis{T1,D}() where {T1,D}
        stencils = ntuple(_ -> UnitRange{T1}(-1:2), D)
        NN       = prod(length.(stencils))
        return new{T1,D,NN}(stencils)
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
        error("invalid B-spline node type")   # constant text: a GPU kernel cannot build a string
    end   
    return ϕ,∂ϕ*Δx⁻¹ 
end

@inline function eval_basis(mpts::Point{T1,T2,1,CM}, mesh::Mesh{T1,T2,1}, basis::Basis{T1,T2,1,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,CM}
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
@inline function eval_basis(mpts::Point{T1,T2,2,CM}, mesh::Mesh{T1,T2,2}, basis::Basis{T1,T2,2,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,CM}
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
@inline function eval_basis(mpts::Point{T1,T2,3,CM}, mesh::Mesh{T1,T2,3}, basis::Basis{T1,T2,3,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:BSplineBasis,CM}
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
