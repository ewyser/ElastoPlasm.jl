@kernel inbounds = true function flip_nd_n2p(mp,mesh,dt)
    p = @index(Global)
    if p≤mp.nmp    
        # flip update
        for dim ∈ 1:mesh.dim
            δa = δv = 0.0
            for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
                δa += (mp.ϕ∂ϕ[nn,p,1]*mesh.a[dim,no])
                δv += (mp.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no])
            end
            mp.v[dim,p]+= dt*δa 
            mp.x[dim,p]+= dt*δv
            # find maximum velocity component over mps
            @atom mp.vmax[dim] = max(mp.vmax[dim],abs(mp.v[dim,p]))
        end
    end  
end
@kernel inbounds = true function pic_nd_n2p(mp,mesh,dt)
    p = @index(Global)
    if p≤mp.nmp    
        for dim ∈ 1:mesh.dim
            δv = 0.0
            # pic update
            for (nn,no) ∈ enumerate(mp.p2n[:,p]) if no<1 continue end
                δv += mp.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no]
            end
            mp.v[dim,p] = δv 
            mp.x[dim,p]+= dt*δv
            # find maximum velocity component over mps
            @atom mp.vmax[dim] = max(mp.vmax[dim],abs(mp.v[dim,p]))
        end
    end  
end
function n2p(mp,mesh,dt,instr)
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(ndrange=mp.nmp,mp,mesh,dt);sync(CPU())
    return nothing
end