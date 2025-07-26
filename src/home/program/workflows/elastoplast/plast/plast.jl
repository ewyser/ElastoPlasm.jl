function init_plast(instr)
    kernel1 = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mp,ϵpII,cmp,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        kernel2 = DP!(CPU())
    elseif instr[:plast][:constitutive] == "J2"
        kernel2 = J2!(CPU())
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mp,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(cmpr[:cmType])"))
    end 
    return (;nonloc! = kernel1, retmap! = kernel2,) 
end
function plast(mp,mesh,cmpr,instr)
    if instr[:plast][:status] 
        # nonlocal regularization
        if instr[:nonloc][:status]
            ls      = instr[:nonloc][:ls]
            mp.e2p.= Int(0)
            mp.p2p.= Int(0)
            mp.ϵpII[2,:].= 0.0
            W,w     = spzeros(mp.nmp),spzeros(mp.nmp,mp.nmp)
            for proc ∈ ["tplgy","p->q","p<-q"]
                instr[:cairn][:elastoplast][:plast].nonloc!(ndrange=mp.nmp,W,w,mp,mesh,ls,proc);sync(CPU())
            end
        end
        # plastic return-mapping dispatcher
        instr[:cairn][:elastoplast][:plast].retmap!(ndrange=mp.nmp,mp,cmpr,instr);sync(CPU())
        ηmax = 0
    else 
        ηmax = 0 
    end

    return ηmax::Int64
end