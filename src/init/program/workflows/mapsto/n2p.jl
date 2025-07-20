@kernel inbounds = true function flip_nd_n2p(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp    
        # flip update
        for dim ∈ 1:meD.nD
            δa = δv = 0.0
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                δa += (mpD.ϕ∂ϕ[nn,p,1]*meD.an[no,dim])
                δv += (mpD.ϕ∂ϕ[nn,p,1]*meD.vn[no,dim])
            end
            mpD.v[p,dim]+= Δt*δa 
            mpD.x[p,dim]+= Δt*δv
        end
    end  
end
@kernel inbounds = true function pic_nd_n2p(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp    
        for dim ∈ 1:meD.nD
            δv = 0.0
            # pic update
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                δv += mpD.ϕ∂ϕ[nn,p,1]*meD.vn[no,dim]
            end
            mpD.v[p,dim] = δv 
            mpD.x[p,dim]+= Δt*δv
        end
    end  
end
function n2p(mpD,meD,Δt,instr)
    # mapping to material point
    instr[:cairn][:mapsto][:map].n2p!(ndrange=mpD.nmp,mpD,meD,Δt);sync(CPU())
    return nothing
end