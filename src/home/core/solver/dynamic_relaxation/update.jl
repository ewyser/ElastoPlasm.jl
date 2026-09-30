"""
    update(mpts, mesh, basis)

End of a converged load step: accept `F`, the logarithmic strain and the stress as the new
converged state, and move each material point by `Δu_p = Σ S_vp Δu_v` (Sutcliffe et al.
2026, Eq 23). Volume, porosity and density follow `deform.jl`, so the solid mass `ρ·Ω`
stays equal to `mpts.s.m`.
"""
@kernel inbounds = true function update(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        mpts.s.Fn[p] = mpts.s.Fᵢⱼ[p]
        mpts.s.ϵn[p] = mpts.s.ϵᵢ[p]
        Δu = zero(SVector{D,T2})
        for nn ∈ 1:mesh.prprt.nn
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            Δu += basis.N[nn,p] * mesh.s.v[no]
        end
        mpts.s.u[p] = Δu
        mpts.x[p]  += Δu
        J             = det(mpts.s.Fᵢⱼ[p])
        mpts.J[p]     = J
        mpts.s.σᵢ[p] = CauchyStress(mpts.s.τᵢ[p], T2(1.0)/J)
        mpts.Ω[p]     = J * mpts.Ω₀[p]
        n             = T2(1.0) - T2(1.0)/J*(T2(1.0)-mpts.n₀[p])
        mpts.s.ρ[p]   = (T2(1.0)-n)*mpts.s.ρ₀[p]
        mpts.n[p]     = min(max(n, T2(0.0)), T2(1.0))
    end
end
