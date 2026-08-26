@kernel inbounds = true function update(mpts::Point{T1,T2,2,CM,TM,TV,TS,ST},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},dt::T2) where {T1,T2,CM,TM,TV,TS,ST<:InfinitesimalStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        # update material point displacement
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no, N, _ = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            δvx += N * mesh.s.v[no][1]
            δvy += N * mesh.s.v[no][2]
        end
        mpts.J[p] = det(mpts.s.Fᵢⱼ[p])
        mpts.Ω[p] = mpts.J[p] * mpts.Ω₀[p]
        mpts.s.ρ[p] = mpts.s.ρ[p] / mpts.ΔJ[p]
        mpts.n[p]   = T2(1.0) - T2(1.0) / mpts.J[p] * (T2(1.0) - mpts.n[p])
        mpts.s.u[p] = SVector{2,T2}(dt * δvx, dt * δvy)
        mpts.x[p] += mpts.s.u[p]
        mpts.s.Fn[p] = mpts.s.Fᵢⱼ[p]
        mpts.s.σn[p] = mpts.s.σᵢⱼ[p]
    end
end

@kernel inbounds = true function update(mpts::Point{T1,T2,2,CM,TM,TV,TS,ST},mesh::Mesh{T1,T2,2},basis::Basis{T1,T2,2},dt::T2) where {T1,T2,CM,TM,TV,TS,ST<:LogarithmicStrain}
    p = @index(Global)
    if p ≤ mpts.nmp
        # store converged current deformation gradient and logarithmic strain
        mpts.s.Fn[p] = mpts.s.Fᵢⱼ[p]
        mpts.s.ϵn[p] = mpts.s.ϵᵢⱼ[p]
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            no, N, _ = eval_basis(mpts, mesh, basis, p, nn)
            if iszero(no) continue end
            δvx += N * mesh.s.v[no][1]
            δvy += N * mesh.s.v[no][2]
        end
        mpts.s.u[p] = SVector{2,T2}(dt * δvx, dt * δvy)
        mpts.x[p] += mpts.s.u[p]
        mpts.J[p] = det(mpts.s.Fn[p])
        mpts.s.σᵢⱼ[p] = CauchyStress(mpts.s.τᵢⱼ[p], T2(1.0)/mpts.J[p])

        # update material point's volume        
        mpts.Ω[p] = mpts.J[p] * mpts.Ω₀[p]

        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p] = mpts.s.ρ₀[p] / mpts.J[p] 
        mpts.n[p]   = T2(1.0) - T2(1.0) / mpts.J[p] * (T2(1.0) - mpts.n₀[p])
    end
end