# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Concrete MLSBasis type — Moving Least Squares shape functions (structured grid)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# See Cao et al. 2025, "Unstructured moving least squares material point methods:
# a stable kernel approach with continuous gradient reconstruction on general
# unstructured tessellations", Comput. Mech. 75:655-678, §2.3.1 (structured-grid
# MLS-MPM). Only the structured-grid case is implemented — ElastoPlasm's Mesh is
# always structured, so the paper's unstructured diminishing-function machinery
# (§2.3.2-2.3.5) does not apply here.
#
# At a particle p, with neighbor nodes x_i and raw (uncorrected) kernel weights w_i:
#   p_i  = [1; x_i - x_p]                      (linear polynomial basis, D+1 vector)
#   M_p  = Σ_i w_i p_i p_iᵀ                     ((D+1)x(D+1) moment matrix)
#   N_i  = w_i (M_p⁻¹ p_i)[1]                   (shape value)
#   ∂N_i = w_i (M_p⁻¹ p_i)[2:end]                (spatial gradient ∇_z w_i(z;x_p)|_{z=x_p})
#
# M_p is built fresh inside eval_basis from all NN neighbors of the particle — cheap
# because eval_basis is only called from the once-per-particle-per-timestep `shpfun!`
# ignite kernel (see ignite/misc.jl), not from every P2G/G2P/update call site.
#
# Unlike BSplineBasis's t1_ϕ∂ϕ..t4_ϕ∂ϕ (node-type-dependent boundary correction), the
# raw kernel here is a single universal formula — MLS's moment-matrix projection is
# what restores partition-of-unity/consistency near boundaries, so no node-type
# classification (basis.type) is needed for this basis kind.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export MLSBasis

struct MLSBasis{T,D,NN} <: AbstractBasis
    stencils::NTuple{D,UnitRange{T}}
    #
    function MLSBasis{T,D}() where {T,D}
        stencils = ntuple(_ -> UnitRange{T}(-1:2), D)
        NN       = prod(length.(stencils))
        return new{T,D,NN}(stencils)
    end
end

"""Universal quadratic B-spline kernel value, support radius 1.5h. ξ = δx/h (either sign)."""
@inline function bspline(ξ::T2) where {T2}
    a = abs(ξ)
    if a < T2(0.5)
        return T2(0.75) - a^2
    elseif a < T2(1.5)
        return T2(0.5) * (T2(1.5) - a)^2
    else
        return T2(0.0)
    end
end

@inline function eval_basis(mpts::Point{T1,T2,D,E,R}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D,NN,K}, ip::T1, nn::T1) where {T1,T2,D,NN,K<:MLSBasis,E,R}
    p2n = basis.p2n[ip]
    no  = p2n[nn]
    if iszero(no)
        return no, T2(0.0), zero(SVector{D,T2})
    end
    h  = mesh.prprt.h
    xp = mpts.x[ip]
    D1 = D + 1
    # raw kernel weights and polynomial basis [1;δx] for every candidate neighbor
    w = zero(MVector{NN,T2})
    P = ntuple(NN) do k
        nk = p2n[k]
        if iszero(nk)
            return zero(SVector{D1,T2})
        end
        δx = mesh.x[nk] - xp
        wk = one(T2)
        for d ∈ 1:D
            wk *= bspline(δx[d]/h[d]) / h[d]
        end
        w[k] = wk
        SVector{D1,T2}(one(T2), δx...)
    end
    # moment matrix and its inverse
    M = zero(MMatrix{D1,D1,T2})
    for k ∈ 1:NN
        M .+= w[k] .* (P[k] * P[k]')
    end
    Minv = inv(SMatrix(M))
    # shape value/gradient for the requested neighbor nn
    Mp = Minv * P[nn]
    N  = w[nn] * Mp[1]
    ∂N = w[nn] * SVector{D,T2}(ntuple(d -> Mp[d+1], D))
    return no, N, ∂N
end
