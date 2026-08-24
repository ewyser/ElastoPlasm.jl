# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete LinearBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export LinearBasis

struct LinearBasis{T,D,NN} <: AbstractBasis
    stencils::NTuple{D,UnitRange{T}}
    #
    function LinearBasis{T,D}() where {T,D}
        stencils = ntuple(_ -> UnitRange{T}(0:1), D)
        NN       = prod(length.(stencils))
        return new{T,D,NN}(stencils)
    end
end

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


@inline function eval_basis(mpts::Point{T1,T2,1,E}, mesh::Mesh{T1,T2,1}, basis::Basis{T1,T2,1,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:LinearBasis,E}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{1,T2}(0.0)
    else
        ϕξ,∂ϕξ = N∂N((mpts.x[ip][1]-mesh.x[no][1]),mesh.prprt.h[1])
        # convolution of basis function
        N, ∂N  = T2(ϕξ), SVector{1,T2}(∂ϕξ)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,2,E}, mesh::Mesh{T1,T2,2}, basis::Basis{T1,T2,2,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:LinearBasis,E}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{2,T2}(0.0, 0.0)
    else
        ϕξ,∂ϕξ = N∂N((mpts.x[ip][1]-mesh.x[no][1]),mesh.prprt.h[1])
        ϕη,∂ϕη = N∂N((mpts.x[ip][2]-mesh.x[no][2]),mesh.prprt.h[2])
        # convolution of basis function
        N, ∂N  = T2(ϕξ*ϕη), SVector{2,T2}(∂ϕξ*ϕη, ϕξ*∂ϕη)
    end
    return no, N, ∂N
end
@inline function eval_basis(mpts::Point{T1,T2,3,E}, mesh::Mesh{T1,T2,3}, basis::Basis{T1,T2,3,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:LinearBasis,E}
    no = basis.p2n[ip][nn]
    if iszero(no)
        N, ∂N = T2(0.0), SVector{3,T2}(0.0, 0.0, 0.0)
    else
        ϕξ,∂ϕξ = N∂N((mpts.x[ip][1]-mesh.x[no][1]),mesh.prprt.h[1])
        ϕη,∂ϕη = N∂N((mpts.x[ip][2]-mesh.x[no][2]),mesh.prprt.h[2])
        ϕζ,∂ϕζ = N∂N((mpts.x[ip][3]-mesh.x[no][3]),mesh.prprt.h[3])
        # convolution of basis function
        N, ∂N  = T2(ϕξ*ϕη*ϕζ), SVector{3,T2}(∂ϕξ*ϕη*ϕζ, ϕξ*∂ϕη*ϕζ, ϕξ*ϕη*∂ϕζ)
    end
    return no, N, ∂N
end