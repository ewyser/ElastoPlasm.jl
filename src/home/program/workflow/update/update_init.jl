function init_update(instr::NamedTuple; cairn::NamedTuple=(;))
    if instr[:fwrk][:deform] == "finite"
        cairn = merge(cairn, (; deform! = finite_deform(CPU()), ))
    elseif instr[:fwrk][:deform] == "infinitesimal"
        cairn = merge(cairn, (; deform! = infinitesimal_deform(CPU()), ))
    end
    cairn = merge(cairn, (; heat! = heatflux(CPU()), ))
    if instr[:basis][:which] == "gimpm"
        cairn = merge(cairn, (; domain! = undeformed(CPU()), ))
        if instr[:fwrk][:deform] == "finite"
            if instr[:basis][:how] == "detFij"
                cairn = merge(cairn, (; domain! = detFᵢᵢ(CPU()), ))
            elseif instr[:basis][:how] == "Fii"
                cairn = merge(cairn, (; domain! = Fᵢᵢ(CPU()), ))
            elseif instr[:basis][:how] == "Uii"
                cairn = merge(cairn, (; domain! = Uᵢᵢ(CPU()), ))
            end
        elseif instr[:fwrk][:deform] == "infinitesimal"
            if instr[:basis][:how] == "detΔFij"
                cairn = merge(cairn, (; domain! = detΔFᵢᵢ(CPU()), ))
            elseif instr[:basis][:how] == "ΔFii"
                cairn = merge(cairn, (; domain! = ΔFᵢᵢ(CPU()), ))
            elseif instr[:basis][:how] == "ΔUii"
                cairn = merge(cairn, (; domain! = ΔUᵢᵢ(CPU()), ))
            end
        end
    end
    if instr[:fwrk][:locking]
        cairn = merge(cairn, (; ΔJn! = ΔJn(CPU()), ΔJs! = ΔJs(CPU()), ΔJp! = ΔJp(CPU()),))
    end
    if !instr[:perf][:status]
        if instr[:fwrk][:deform] == "finite"
            cairn = merge(cairn, (; elast! = finite_elast(CPU()), ))
        elseif instr[:fwrk][:deform] == "infinitesimal"
            cairn = merge(cairn, (; elast! = infinitesimal_elast(CPU()), ))
        end
    end
    cairn = merge(cairn, (; nonloc! = nonlocal(CPU()), ))
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mpts,ϵpII,cmp,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        if instr[:fwrk][:deform] == "finite"
            cairn = merge(cairn, (; retmap! = finite_DP(CPU()), ))
        elseif instr[:fwrk][:deform] == "infinitesimal"
            cairn = merge(cairn, (; retmap! = infinitesimal_DP(CPU()), ))
        end
    elseif instr[:plast][:constitutive] == "J2"
        if instr[:fwrk][:deform] == "finite"
            cairn = merge(cairn, (; retmap! = finite_J2(CPU()), ))
        elseif instr[:fwrk][:deform] == "infinitesimal"
            cairn = merge(cairn, (; retmap! = infinitesimal_J2(CPU()), ))
        end
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mpts,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(instr[:plast][:constitutive])"))
    end
    
    return cairn
end