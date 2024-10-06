@kernel inbounds = true function kernel_momentum(mpD,meD)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for nn ∈ 1:meD.nn
                @atom meD.pn[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
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
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            mpD.u[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end
    end
end
function DM!(mpD,meD,Δt)
    # initialize for DM
    meD.pn.= 0.0
    meD.vn.= 0.0
    # accumulate material point contributions
    @isdefined(DMp2n!) ? nothing : DMp2n! = kernel_momentum(CPU())
    DMp2n!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
    # solve for nodal incremental displacement
    @isdefined(DMsolve!) ? nothing : DMsolve! = kernel_velocity(CPU())
    DMsolve!(meD; ndrange=meD.nno[end]);sync(CPU())
    # update material point's displacement
    @isdefined(DMdispl!) ? nothing : DMdispl! = kernel_displacement(CPU())
    DMdispl!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
    return nothing
end