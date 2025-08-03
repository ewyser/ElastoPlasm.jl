@kernel inbounds = true function flip_nd_n2p(mp::Point{T1,T2},mesh,dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp    
        # flip update
        for dim ∈ 1:mesh.dim
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if iszero(no) continue end
                mp.s.v[dim,p]+= dt*(mp.ϕ∂ϕ[nn,p,1]*mesh.a[dim,no])
                mp.x[dim,p]  += dt*(mp.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no])
            end
            # find maximum velocity component over mps
            @atom mp.vmax[dim] = max(mp.vmax[dim],abs(mp.s.v[dim,p]))
        end
    end  
end
@kernel inbounds = true function pic_nd_n2p(mp::Point{T1,T2},mesh,dt::T2) where {T1,T2}
    p = @index(Global)
    if p≤mp.nmp    
        # pic update
        for dim ∈ 1:mesh.dim
            δv = 0.0
            for nn ∈ 1:mesh.nn
                no = mp.p2n[nn,p]
                if iszero(no) continue end
                δv += mp.ϕ∂ϕ[nn,p,1]*mesh.v[dim,no]
            end
            mp.s.v[dim,p] = δv 
            mp.x[dim,p]  += dt*δv
            # find maximum velocity component over mps
            @atom mp.vmax[dim] = max(mp.vmax[dim],abs(mp.s.v[dim,p]))
        end
    end  
end
function n2p(mp::Point{T1,T2},mesh,dt::T2,instr::Dict) where {T1,T2}
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(ndrange=mp.nmp,mp,mesh,dt);sync(CPU())
    return nothing
end