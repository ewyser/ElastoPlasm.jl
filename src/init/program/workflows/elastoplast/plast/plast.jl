function init_plast(instr)
    kernel1 = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mpD,ϵpII,cmParam,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        kernel2 = DP!(CPU())
    elseif instr[:plast][:constitutive] == "J2"
        kernel2 = J2!(CPU())
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mpD,cmParam,instr[:fwrk])
    else
        err_msg = "$(cmParam[:cmType]): invalid return mapping for plastic correction"
        throw(error(err_msg))
    end 
    return (;nonloc! = kernel1, retmap! = kernel2,) 
end
function plast!(mpD,meD,cmParam,instr)
    # nonlocal regularization
    if cmParam[:nonlocal][:status]
        ls      = cmParam[:nonlocal][:ls]
        mpD.e2p.= Int(0)
        mpD.p2p.= Int(0)
        mpD.ϵpII[:,2].= 0.0
        W,w     = spzeros(mpD.nmp),spzeros(mpD.nmp,mpD.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr[:cairn][:elastoplast][:plast].nonloc!(W,w,mpD,meD,ls,proc; ndrange=mpD.nmp);sync(CPU())
        end
    end
    # plastic return-mapping dispatcher
    instr[:cairn][:elastoplast][:plast].retmap!(mpD,cmParam,instr; ndrange=mpD.nmp);sync(CPU())
    ηmax = 0
    return ηmax::Int64
end