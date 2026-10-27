# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete LinearBasis types
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export GimpBasis

# how GIMP updates a particle's domain (`basis.domain`, `particle_domain` in `common/domain.jl`)
abstract type AbstractDomain end
"The stretch update, `ℓ = ℓ₀ ∘ diag(U)` with `U = √(FᵀF)` (Charlton et al. 2017; AMPLE, Coombs & Augarde 2020 Eq 34)."
struct StretchDomain <: AbstractDomain end
"uGIMP: the domain keeps its initial size (Bardenhagen & Kober 2004)."
struct UndeformedDomain <: AbstractDomain end
"""
    get_domain(which::String) -> AbstractDomain

Map a solver's `basis.domain` string (`"stretch"`, `"undeformed"`) to its domain update.
"""
function get_domain(which::String)
    which == "stretch"    && return StretchDomain()
    which == "undeformed" && return UndeformedDomain()
    throw(ArgumentError("basis.domain = \"$which\" is not one of \"stretch\", \"undeformed\""))
end

struct GimpBasis{T,D,NN,P<:AbstractDomain} <: AbstractBasis
    stencils::NTuple{D,UnitRange{T}}
    domain  ::P                                              # how a particle's domain is updated
    #
    function GimpBasis{T,D}(domain::P) where {T,D,P<:AbstractDomain}
        stencils = ntuple(_ -> UnitRange{T}(-1:2), D)
        NN       = prod(length.(stencils))
        return new{T,D,NN,P}(stencils, domain)
    end
end

@inline function S∂S(δx::T2,h::T2,lp::T2) where {T2}
    # Follows AMPLE-m's SvpGIMP.m (Bardenhagen & Kober 2004); region labels A-E as in the AMPLE paper.
    # δx is the signed (particle - node) offset, matching AMPLE's `xp-xv`.
    if (-h-lp) < δx <= (-h+lp)                    # A: partial overlap of domain and element
        S  = (h+lp+δx)^2/(T2(4.0)*h*lp)
        ∂S = (h+lp+δx)/(T2(2.0)*h*lp)
    elseif (-h+lp) < δx <= -lp                    # B: full overlap of domain and element (left)
        S  = T2(1.0)+δx/h
        ∂S = T2(1.0)/h
    elseif -lp < δx <= lp                         # C: partial overlap of domain and two elements
        S  = T2(1.0)-(δx^2+lp^2)/(T2(2.0)*h*lp)
        ∂S = -δx/(h*lp)
    elseif lp < δx <= (h-lp)                      # D: full overlap of domain and element (right)
        S  = T2(1.0)-δx/h
        ∂S = -T2(1.0)/h
    elseif (h-lp) < δx <= (h+lp)                  # E: partial overlap of domain and element
        S  = (h+lp-δx)^2/(T2(4.0)*h*lp)
        ∂S = -(h+lp-δx)/(T2(2.0)*h*lp)
    else                                          # zero overlap
        S,∂S = T2(0.0),T2(0.0)
    end
    return S,∂S
end
@inline function eval_basis(mpts::Point{T1,T2,1,CM}, mesh::Mesh{T1,T2,1}, basis::Basis{T1,T2,1,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,CM}
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
@inline function eval_basis(mpts::Point{T1,T2,2,CM}, mesh::Mesh{T1,T2,2}, basis::Basis{T1,T2,2,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,CM}
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
@inline function eval_basis(mpts::Point{T1,T2,3,CM}, mesh::Mesh{T1,T2,3}, basis::Basis{T1,T2,3,NN,K}, ip::T1, nn::T1) where {T1,T2,NN,K<:GimpBasis,CM}
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