@kernel inbounds = true function update(mpts::Point{T1,T2,D,B,E,R},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D<:TwoDimension,B<:AbstractBasis,E<:LinearElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        # update material point displacement
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N * mesh.s.v[1,no]
            δvy += N * mesh.s.v[2,no]
        end
        mpts.s.u[1,p] = dt * δvx
        mpts.s.u[2,p] = dt * δvy

        # update material point's coordinates
        mpts.x[1,p] += dt * δvx
        mpts.x[2,p] += dt * δvy

        # update material point's deformation jacobian
        mpts.J[p] = det(mpts.s.Fᵢⱼ[:,:,p])

        # update material point's volume        
        mpts.Ω[p] = mpts.J[p] * mpts.Ω₀[p]

        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p] = mpts.s.ρ[p] / mpts.ΔJ[p]     
        mpts.n[p]   = T2(1.0) - T2(1.0) / mpts.J[p] * (T2(1.0) - mpts.n[p])

        # update deformation gradient and cauchy stress
        mpts.s.Fn[:,:,p] = copy(mpts.s.Fᵢⱼ[:,:,p])
        mpts.s.σn[:,p]   = copy(mpts.s.σᵢ[:,p])
    end
end

@kernel inbounds = true function update(mpts::Point{T1,T2,D,B,E,R},mesh::Mesh{T1,T2,D},dt::T2) where {T1,T2,D<:TwoDimension,B<:AbstractBasis,E<:FiniteElasticity,R<:AbstractRheology}
    p = @index(Global)
    if p ≤ mpts.nmp
        # store converged current deformation gradient and logarithmic strain
        mpts.s.Fn[:,:,p] = mpts.s.Fᵢⱼ[:,:,p]
        mpts.s.ϵn[:,:,p] = mpts.s.ϵᵢⱼ[:,:,p]

        # update material point displacement
        δvx = δvy = T2(0.0)
        for nn ∈ 1:mesh.prprt.nn
            # buffering & compute basis functions on-the-fly
            no, N, _ = basis(mpts, mesh, p, nn)
            if iszero(no) continue end
            δvx += N * mesh.s.v[1,no]
            δvy += N * mesh.s.v[2,no]
        end
        mpts.s.u[1,p] = dt * δvx
        mpts.s.u[2,p] = dt * δvy

        # update material point's coordinates
        mpts.x[1,p] += dt * δvx
        mpts.x[2,p] += dt * δvy

        # update material point's deformation jacobian
        mpts.J[p] = det(mpts.s.Fn[:,:,p])

        # update material point's cauchy stress
        mpts.s.σᵢ[:,p] = mpts.s.τᵢ[:,p]./mpts.J[p]

        # update material point's volume        
        mpts.Ω[p] = mpts.J[p] * mpts.Ω₀[p]

        # update material point density and positivity-preserving porosity
        mpts.s.ρ[p] = mpts.s.ρ₀[p] / mpts.J[p] 
        mpts.n[p]   = T2(1.0) - T2(1.0) / mpts.J[p] * (T2(1.0) - mpts.n₀[p])
    end
end