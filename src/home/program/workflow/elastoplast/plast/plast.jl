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
function plast(mp,mesh,cmpr,instr)
    # nonlocal regularization
    if instr[:nonloc][:status]
        ls             = instr[:nonloc][:ls]
        mp.e2p        .= Int(0)
        mp.p2p        .= Int(0)
        mp.s.ϵpII[2,:].= 0.0
        W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr[:cairn][:elastoplast][:plast].nonloc!(ndrange=mp.nmp,W,w,mp,mesh,ls,proc);sync(CPU())
        end
    else
        mp.s.ϵpII[2,:].= mp.s.ϵpII[1,:]
    end
    # plastic return-mapping dispatcher
    instr[:cairn][:elastoplast][:plast].retmap!(ndrange=mp.nmp,mp,cmpr);sync(CPU())
    ηmax = 0
    return ηmax::Int64
end