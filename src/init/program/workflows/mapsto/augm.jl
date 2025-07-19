@kernel inbounds = true function augm_momentum(mpD,meD)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
            end
        end
    end
end
@kernel inbounds = true function augm_velocity(meD)
    n = @index(Global)
    for dim ∈ 1:meD.nD
        if n≤meD.nno[end] 
            if meD.mn[n]>0.0
                meD.vn[n,dim] = (meD.pn[n,dim]*(1.0/meD.mn[n])*meD.bc[n,dim])
            end   
        end
    end
end
@views @kernel inbounds = true function augm_displacement(mpD,meD,Δt)
    p = @index(Global)
    # flip update
    for dim ∈ 1:meD.nD
        Δu = 0.0
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            Δu += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[no,dim])
        end
        mpD.u[p,dim]+= Δu
    end

end