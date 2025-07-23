function init_plast(instr)
    kernel1 = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mp,ϵpII,cmParam,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        kernel2 = DP!(CPU())
    elseif instr[:plast][:constitutive] == "J2"
        kernel2 = J2!(CPU())
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mp,cmParam,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(cmParam[:cmType])"))
    end 
    return (;nonloc! = kernel1, retmap! = kernel2,) 
end
function plast(mp,mesh,cmParam,instr)
    if instr[:plast][:status] 
        # nonlocal regularization
        if cmParam[:nonlocal][:status]
            ls      = cmParam[:nonlocal][:ls]
            mp.e2p.= Int(0)
            mp.p2p.= Int(0)
            mp.ϵpII[2,:].= 0.0
            W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
            for proc ∈ ["tplgy","p->q","p<-q"]
                instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mp,mesh,ls,proc; ndrange=mp.nmp);sync(CPU())
            end
        end
        # plastic return-mapping dispatcher
        instr[:cairn][:elastoplast][:plast].retmap!(mp,cmParam,instr; ndrange=mp.nmp);sync(CPU())
        ηmax = 0
    else 
        ηmax = 0 
    end

    return ηmax::Int64
end