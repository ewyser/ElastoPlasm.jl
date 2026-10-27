"""
    dr_kinematics(mpts, mesh, basis)

Trial kinematics of one dynamic-relaxation iteration (Sutcliffe et al. 2026, Eqs 14-15):
`ΔF = I + Σ Δuᵥ ⊗ ∂Sᵥₚ/∂X̂`, `F = ΔF·Fₙ`, with the nodal displacement increment `Δu`
held in `mesh.s.v` and `∂S/∂X̂` the shape-function gradients cached by `ignite` at the
start of the load step. Resets `ϵᵢ ← ϵn`, so that `elast!` pushes the *converged* strain
forward (Eq 16) and no state carries over from one iteration to the next.
"""
@kernel inbounds = true function dr_kinematics(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        TM = eltype(mpts.s.ΔFᵢⱼ)
        ∇u = zeros(MMatrix{D,D,T2})
        for nn ∈ 1:stencil_size(basis)
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            ∂N = ∂Nrow(basis.∂N, nn, p, Val(D))
            Δu = mesh.s.v[no]
            for i ∈ 1:D, j ∈ 1:D
                ∇u[i,j] += Δu[i]*∂N[j]
            end
        end
        ΔFᵢⱼ           = TM(I) + TM(∇u)
        mpts.s.ΔFᵢⱼ[p] = ΔFᵢⱼ
        mpts.s.Fᵢⱼ[p]  = ΔFᵢⱼ * mpts.s.Fn[p]
        mpts.s.ϵᵢ[p]  = mpts.s.ϵn[p]
    end
end

"""
    dr_fint(mpts, mesh, basis, g)

Accumulate the out-of-balance force `fᵉˣᵗ − fⁱⁿᵗ` into `mesh.s.oobf` (Sutcliffe et al. 2026,
Eqs 19-22): `fⁱⁿᵗ = Σ ∇ₓS·τ/det(F)·V_p` with `∂S/∂x = ∂S/∂X̂·ΔF⁻¹` (Eq 21) and
`V_p = Vₙ·det(ΔF)` (Eq 20), and `fᵉˣᵗ = Σ S·m_p·g`.
"""
@kernel inbounds = true function dr_fint(mpts::Point{T1,T2,D},mesh::Mesh{T1,T2,D},basis::Basis{T1,T2,D},g::Vector{T2}) where {T1,T2,D}
    p = @index(Global)
    if p ≤ mpts.nmp
        ΔFᵢⱼ = mpts.s.ΔFᵢⱼ[p]
        ΔF⁻¹ = inv(ΔFᵢⱼ)
        σᵢⱼ  = get_tensor(mpts.s.τᵢ[p], Val(D)) / det(mpts.s.Fᵢⱼ[p])
        Ω    = mpts.Ω[p] * det(ΔFᵢⱼ)
        m    = mpts.s.m[p]
        for nn ∈ 1:stencil_size(basis)
            no = basis.p2n[p][nn]
            if iszero(no) continue end
            N  = basis.N[nn,p]
            ∂N = ΔF⁻¹' * ∂Nrow(basis.∂N, nn, p, Val(D))
            f  = Ω * (σᵢⱼ * ∂N)
            for d ∈ 1:D
                @atom mesh.s.oobf[d,no] += N * m * g[d] - f[d]
            end
        end
    end
end
