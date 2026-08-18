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

"""Raw kernel weights, polynomial basis [1;δx], and inverted moment matrix for particle ip —
the once-per-particle work shared by every neighbor's shape value/gradient."""
@inline function mls_moments(mpts::Point{T1,T2,D,E,R}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D,NN,K}, ip::T1) where {T1,T2,D,NN,K<:MLSBasis,E,R}
    p2n = basis.p2n[ip]
    h   = mesh.prprt.h
    xp  = mpts.x[ip]
    D1  = D + 1
    # weight + polynomial basis [1;δx] per neighbor — plain scalar writes into fixed-size
    # mutable buffers (no closures) to avoid StaticArrays/ntuple-do-block allocation pitfalls
    w = zero(MVector{NN,T2})
    P = zero(MMatrix{NN,D1,T2})
    for k ∈ 1:NN
        nk = p2n[k]
        if iszero(nk) continue end
        δx = mesh.x[nk] - xp
        wk = one(T2)
        for d ∈ 1:D
            wk *= bspline(δx[d]/h[d]) / h[d]
        end
        w[k]   = wk
        P[k,1] = one(T2)
        for d ∈ 1:D
            P[k,d+1] = δx[d]
        end
    end
    M = zero(MMatrix{D1,D1,T2})
    for k ∈ 1:NN
        pk = SVector{D1,T2}(ntuple(d -> P[k,d], Val(D1)))
        M .+= w[k] .* (pk * pk')
    end
    return SVector{NN,T2}(w), P, inv(SMatrix(M))
end

@inline function eval_basis(mpts::Point{T1,T2,D,E,R}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D,NN,K}, ip::T1, nn::T1) where {T1,T2,D,NN,K<:MLSBasis,E,R}
    no = basis.p2n[ip][nn]
    if iszero(no)
        return no, T2(0.0), zero(SVector{D,T2})
    end
    D1 = D + 1
    w, P, Minv = mls_moments(mpts, mesh, basis, ip)
    pnn = SVector{D1,T2}(ntuple(d -> P[nn,d], Val(D1)))
    Mp  = Minv * pnn
    N   = w[nn] * Mp[1]
    ∂N  = w[nn] * SVector{D,T2}(ntuple(d -> Mp[d+1], Val(D)))
    return no, N, ∂N
end

"""Shape values/gradients for ALL NN neighbors of particle ip in one pass — builds the
(D+1)x(D+1) moment matrix once instead of once per neighbor (used by the MLS-specialized
`shpfun!` kernel, an O(NN) alternative to calling `eval_basis` NN times, O(NN²))."""
@inline function eval_basis_all(mpts::Point{T1,T2,D,E,R}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D,NN,K}, ip::T1) where {T1,T2,D,NN,K<:MLSBasis,E,R}
    D1 = D + 1
    w, P, Minv = mls_moments(mpts, mesh, basis, ip)
    Ns  = zero(MVector{NN,T2})
    ∂Ns = zero(MMatrix{NN,D,T2})
    for k ∈ 1:NN
        pk = SVector{D1,T2}(ntuple(d -> P[k,d], Val(D1)))
        Mp = Minv * pk
        Ns[k] = w[k] * Mp[1]
        for d ∈ 1:D
            ∂Ns[k,d] = w[k] * Mp[d+1]
        end
    end
    return SVector{NN,T2}(Ns), SMatrix{NN,D,T2}(∂Ns)
end

# Overrides the generic per-nn shpfun! (ignite/misc.jl) with the O(NN) batched path above —
# same role (cache basis.N/basis.∂N once per particle per timestep), cheaper for MLS since
# the naive per-nn eval_basis loop would rebuild the moment matrix NN times (O(NN²)).
@kernel inbounds = true function shpfun!(mpts::Point{T1,T2,D,E,R}, mesh::Mesh{T1,T2,D}, basis::Basis{T1,T2,D,NN,K}) where {T1,T2,D,NN,K<:MLSBasis,E,R}
    p = @index(Global)
    if p ≤ mpts.nmp
        N, ∂N = eval_basis_all(mpts, mesh, basis, p)
        basis.N[p]  = N
        basis.∂N[p] = ∂N
    end
end
