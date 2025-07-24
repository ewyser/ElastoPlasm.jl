@kernel inbounds = true function augm_momentum(mp,mesh)
    p = @index(Global)
    for dim ∈ 1:mesh.dim
        if p≤mp.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
                @atom mesh.pn[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.m[p]*mp.v[dim,p])
            end
        end
    end
end
@kernel inbounds = true function augm_velocity(mesh)
    no = @index(Global)
    for dim ∈ 1:mesh.dim
        if no≤mesh.nno[end] 
            if mesh.mn[no]>0.0
                mesh.vn[dim,no] = (mesh.pn[dim,no]*(1.0/mesh.mn[no])*mesh.bc[dim,no])
            end   
        end
    end
end
@views @kernel inbounds = true function augm_displacement(mp,mesh,dt)
    p = @index(Global)
    # flip update
    for dim ∈ 1:mesh.dim
        Δu = 0.0
        for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
            Δu += dt*(mp.ϕ∂ϕ[nn,p,1]*mesh.vn[dim,no])
        end
        mp.u[dim,p]+= Δu
    end
end
function augm(mp,mesh,dt,instr)
    # initialize for DM
    mesh.pn.= 0.0
    mesh.vn.= 0.0
    # accumulate material point contributions
    instr[:cairn][:mapsto][:augm].p2n!(ndrange=mp.nmp,mp,mesh);sync(CPU())
    # solve for nodal incremental displacement
    instr[:cairn][:mapsto][:augm].solve!(ndrange=mesh.nno[end],mesh);sync(CPU())
    # update material point's displacement
    instr[:cairn][:mapsto][:augm].Δu!(ndrange=mp.nmp,mp,mesh,dt);sync(CPU())
    return nothing
end