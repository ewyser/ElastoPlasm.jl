function init_plast(instr)
    kernel1 = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mp,ϵpII,cmp,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        if instr[:fwrk][:deform] == "finite"
            kernel2 = finite_DP(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            kernel2 = infinitesimal_DP(CPU())
        end
    elseif instr[:plast][:constitutive] == "J2"
        if instr[:fwrk][:deform] == "finite"
            kernel2 = finite_J2(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            kernel2 = infinitesimal_J2(CPU())
        end
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mp,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(instr[:plast][:constitutive])"))
    end 
    return (;nonloc! = kernel1, retmap! = kernel2,) 
end
function plast(mp::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,instr::Dict) where {T1,T2}
    # nonlocal regularization
    if instr[:nonloc][:status]
        ls             = T2(instr[:nonloc][:ls])
        mp.e2p        .= T1(0)
        mp.p2p        .= T1(0)
        mp.s.ϵpII[2,:].= T2(0.0)
        W,w     = spzeros(T2,mp.nmp),spzeros(T2,mp.nmp,mp.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr[:cairn][:elastoplast][:plast].nonloc!(ndrange=mp.nmp,W,w,mp,mesh,ls,proc);sync(CPU())
        end
    else
        mp.s.ϵpII[2,:].= mp.s.ϵpII[1,:]
    end
    # plastic return-mapping dispatcher
    instr[:cairn][:elastoplast][:plast].retmap!(ndrange=mp.nmp,mp,cmpr);sync(CPU())
    return nothing
end