@kernel inbounds = true function kernel_momentum(mpD,meD)
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
@kernel inbounds = true function kernel_velocity(meD)
    n = @index(Global)
    for dim ∈ 1:meD.nD
        if n≤meD.nno[end] 
            if meD.mn[n]>0.0
                meD.vn[n,dim] = (meD.pn[n,dim]*(1.0/meD.mn[n])*meD.bc[n,dim])
            end   
        end
    end
end
@views @kernel inbounds = true function kernel_displacement(mpD,meD,Δt)
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
function init_double(trsfr)
    if trsfr=="mUSL"
        p2n!   = kernel_momentum(CPU())
        solve! = kernel_velocity(CPU())
        displ! = kernel_displacement(CPU())
        return (;p2n! = p2n!,solve! = solve!,Δu! = displ!,)
    else
        return nothing
    end    
end
function DM!(mpD,meD,Δt,instr)
    # initialize for DM
    meD.pn.= 0.0
    meD.vn.= 0.0
    # accumulate material point contributions
    instr[:cairn][:augm!].p2n!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
    # solve for nodal incremental displacement
    instr[:cairn][:augm!].solve!(meD; ndrange=meD.nno[end]);sync(CPU())
    # update material point's displacement
    instr[:cairn][:augm!].Δu!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
    return nothing
end