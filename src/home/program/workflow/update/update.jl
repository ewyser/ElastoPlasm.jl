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
    cairn = merge(cairn, (; elast! = elast(CPU()), ))
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

function update(::Val{:elastodynamic!},mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::NamedTuple) where {T1,T2,E,R}
    return elasto(mpts,mesh,cmpr,dt,instr)    
end
function update(::Val{:elastoplastic!},mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::NamedTuple) where {T1,T2,E,R}
    return elastoplast(mpts,mesh,cmpr,dt,instr)    
end
function update(::Val{:thermodynamic!},mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::NamedTuple) where {T1,T2,E,R}
    return thermo(mpts,mesh.t,instr)
end

function elasto(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::Instruction{T1,T2}) where {T1,T2,E,R}
    # update {logarithmic|infinitesimal} strains
    instr.cairn.update.deform!(mpts,mesh.s,dt; ndrange=mpts.nmp);sync(CPU())
    # update material point's domain
    if instr.basis.which == "gimpm"
        instr.cairn.update.domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    if instr.fwrk.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0)/mesh.prprt.dim
        # mapping to mesh 
        instr.cairn.update.ΔJn!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr.cairn.update.ΔJs!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
        # compute determinant Jbar 
        instr.cairn.update.ΔJp!(mpts,mesh,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    instr.cairn.update.elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())
    return nothing
end

function elastoplast(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::Instruction{T1,T2}) where {T1,T2,E,R}
    # update {logarithmic|infinitesimal} strains
    instr.cairn.update.deform!(mpts,mesh.s,dt; ndrange=mpts.nmp);sync(CPU())
    # update material point's domain
    if instr.basis.which == "gimpm"
        instr.cairn.update.domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    if instr.fwrk.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0)/mesh.prprt.dim
        # mapping to mesh 
        instr.cairn.update.ΔJn!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr.cairn.update.ΔJs!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
        # compute determinant Jbar 
        instr.cairn.update.ΔJp!(mpts,mesh,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    instr.cairn.update.elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())
    # plastic corrector
    if instr.nonloc.status
        fill!(mpts.e2p,T1(0))
        fill!(mpts.p2p,T1(0))
        mpts.s.ϵpII[2,:].= T2(0.0)
        W,w     = spzeros(T2,mpts.nmp),spzeros(T2,mpts.nmp,mpts.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr.cairn.update.nonloc!(W,w,mpts,mesh,T2(instr.nonloc.ls),proc; ndrange=mpts.nmp);sync(CPU())
        end
    else
        mpts.s.ϵpII[2,:].= mpts.s.ϵpII[1,:]
    end
    # plastic return-mapping dispatcher
    instr.cairn.update.retmap!(mpts,cmpr; ndrange=mpts.nmp);sync(CPU())
    return nothing
end

function thermo(mpts::Point{T1,T2,E,R},mesh::MeshThermalPhase{T1,T2},instr::Instruction{T1,T2}) where {T1,T2,E,R}
    # update temperature
    instr.cairn.update.heat!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    return nothing
end