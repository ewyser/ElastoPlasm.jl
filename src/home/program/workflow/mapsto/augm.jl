@kernel inbounds = true function augm_momentum(mp::Point{T1,T2},mesh) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            # accumulation
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                @atom mesh.p[dim,no]+= mp.ϕ∂ϕ[nn,p,1]*(mp.s.m[p]*mp.s.v[dim,p])
            end
        end
    end
end
@kernel inbounds = true function augm_velocity(mesh)
    no = @index(Global)
    if no≤mesh.nno[end] 
        if iszero(mesh.mᵢ[no])
            nothing         
        else
            for dim ∈ 1:mesh.dim       
                # apply boundary contidions || forward euler solution
                if mesh.bc[dim,no] == 0.0
                    mesh.v[dim,no] = 0.0                                         
                else
                    mesh.v[dim,no] = mesh.p[dim,no]*(1.0/mesh.mᵢ[no])  
                end
            end
        end
    end
end
@views @kernel inbounds = true function augm_displacement(mp::Point{T1,T2},mesh,dt::T2) where {T1,T2}
    p = @index(Global)
    # flip update
    if p≤mp.nmp
        for dim ∈ 1:mesh.dim 
            Δu = T2(0.0)
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if no < 1 continue end
                Δu += dt*(mp.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no])
            end
            mp.s.u[dim,p]+= Δu
        end
    end
end
function augm(mp::Point{T1,T2},mesh,dt::T2,instr::Dict) where {T1,T2}
    # initialize for DM
    mesh.p.= T2(0.0)
    mesh.v.= T2(0.0)
    # accumulate material point contributions
    instr[:cairn][:mapsto][:augm].p2n!(ndrange=mp.nmp,mp,mesh);sync(CPU())
    # solve for nodal incremental displacement
    instr[:cairn][:mapsto][:augm].solve!(ndrange=mesh.nno[end],mesh);sync(CPU())
    # update material point's displacement
    instr[:cairn][:mapsto][:augm].Δu!(ndrange=mp.nmp,mp,mesh,dt);sync(CPU())
    return nothing
end