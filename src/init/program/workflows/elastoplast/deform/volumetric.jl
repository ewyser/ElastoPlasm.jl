@kernel inbounds = true function ΔJn(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        # accumulation
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            @atom meD.ΔJn[no]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.ΔJ[p])  
        end
    end
end
@kernel inbounds = true function ΔJs(mpD,meD)
    no = @index(Global)
    if no≤meD.nno[end] 
        # solve
        meD.mn[no]>0.0 ? meD.ΔJn[no]/= meD.mn[no] : meD.ΔJn[no] = 0.0
    end
end
@kernel inbounds = true function ΔJp(mpD,meD,dim)
    p = @index(Global)
    if p≤mpD.nmp 
        # mapping back to mp's
        ΔJ = 0.0
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            ΔJ += mpD.ϕ∂ϕ[nn,p,1]*meD.ΔJn[no]/mpD.ΔJ[p]
        end
        @views mpD.ΔFᵢⱼ[:,:,p].*= ΔJ^dim
    end
end
function init_volumetric()
    return (;ΔJn! = ΔJn(CPU()),ΔJs! = ΔJs(CPU()),ΔJp! = ΔJp(CPU()),)
end
function volumetric(mpD,meD,instr)
    if instr[:vollock]
        # init mesh quantities to zero
        meD.ΔJn.= 0.0
        # calculate dimensional cst.
        dim     = 1.0/meD.nD
        # mapping to mesh 
        instr[:cairn][:elastoplast][:Fbar].ΔJn!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr[:cairn][:elastoplast][:Fbar].ΔJs!(mpD,meD; ndrange=meD.nno[end]);sync(CPU())
        # compute determinant Jbar 
        instr[:cairn][:elastoplast][:Fbar].ΔJp!(mpD,meD,dim; ndrange=mpD.nmp);sync(CPU())
    end
    return nothing
end