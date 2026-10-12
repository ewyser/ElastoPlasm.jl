# Nodal solid fraction and exclusion of ill-conditioned nodes (Xie, Navas & López-Querol 2026,
# Engineering with Computers 42:146, §3.4-3.5). A node barely covered by material ("small cut")
# carries a tiny share of mass but a full share of stiffness gradients, which destabilises the
# dynamic relaxation: its slow modes dominate the iteration. `relax!` fixes such nodes (Δu = 0) for
# the load step. Used by the DR solver only: in the explicit solver the fixed nodes feed a spurious
# velocity gradient to nearby particles every step, which accumulates until a particle collapses.
# Fixing also holds back the free surface in proportion to the threshold (see `stab.sfrac`).

"""
    nodal_volume(mpts, mesh, solver) -> Vector

Nodal volume `V_I = Σ_gp S_I,gp V_gp` (Eq 37): the node's share of a completely filled grid,
integrated with quadrature points laid out like the material points (`ni` per direction, taken
from the particles' initial domain) and the problem's own basis functions.
"""
function nodal_volume(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    ni    = round(Int, mesh.prprt.h[1] / (T2(2.0) * mpts.ℓ₀[1][1]))
    mat   = setup_material_constants(solver)
    x     = mpts_populate(mesh.prprt, mat, solver; ni=ni).x
    nmp   = size(x, 2)
    gps   = setup_mpts(mesh, solver, mat; geom=(; xp=x, ni=ni, nmp=nmp, material_fields(mat, nmp)...))
    basis = setup_basis(MechanicalProblem(mesh, gps, setup_time(solver)), solver)
    ignite(gps, mesh, basis, solver)
    V = zeros(T2, mesh.prprt.nno[end])
    for p ∈ 1:gps.nmp, nn ∈ 1:mesh.prprt.nn
        no = basis.p2n[p][nn]
        iszero(no) && continue
        V[no] += basis.N[nn,p] * gps.Ω₀[p]
    end
    return V
end

"""
    solid_fraction(mpts, mesh, basis, r)

Accumulate `Σ_p S_Ip V_p (1 − n_p)/(1 − n₀_p)` into `r`: the solid volume of Eq 36, relative to
each material point's initial solid fraction, so that `r_I / V_I` is 1 inside the material and
falls to 0 across a free surface.
"""
@kernel inbounds = true function solid_fraction(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},r::Vector{T2}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        vs = mpts.Ω[p] * (T2(1.0) - mpts.n[p]) / (T2(1.0) - mpts.n₀[p])
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            @atom r[no] += basis.N[nn,p] * vs
        end
    end
end

"""
    inactive_nodes!(inactive, r, mpts, mesh, basis, V, solver) -> inactive

Mark the nodes that carry material but whose relative solid fraction `r_I / V_I` is below
`solver.stab.sfrac`. `r` is scratch space of the nodal length. Shape functions must be cached
(`ignite`) for the current configuration.
"""
function inactive_nodes!(inactive::Vector{Bool},r::Vector{T2},mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},V::Vector{T2},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    fill!(r, T2(0.0))
    solver.cairn.ignite.sfrac!(mpts,mesh,basis,r; ndrange=mpts.nmp);sync(CPU())
    ε = T2(solver.stab.sfrac)
    @. inactive = (r > T2(0.0)) & (r < ε * V)
    return inactive
end

"""
    node_exclusion(mpts, mesh, solver) -> NamedTuple or nothing

Workspace for excluding ill-conditioned nodes: the nodal volume `V` (Eq 37, computed once for
the fixed background grid) and scratch buffers; `nothing` when `solver.stab.sfrac == 0`.
"""
function node_exclusion(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},solver::S) where {T1,T2,D,S<:AbstractSolver{T1,T2,D}}
    iszero(solver.stab.sfrac) && return nothing
    nno = mesh.prprt.nno[end]
    return (; V=nodal_volume(mpts, mesh, solver), r=zeros(T2, nno), inactive=zeros(Bool, nno))
end
