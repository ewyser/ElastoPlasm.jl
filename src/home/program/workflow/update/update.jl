function init_update(instr::NamedTuple; update::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    if instr[:perf][:status]
        update[:deform!] = deform_fast(CPU())
    else
        update[:deform!] = deform(CPU())
    end
    
    update[:heat!] = heatflux(CPU())
    if instr[:basis][:which] == "gimpm"
        update[:domain!] = undeformed(CPU())
        if instr[:fwrk][:deform] == "finite"
            if instr[:basis][:how] == "detFij"
                update[:domain!] = detFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Fii"
                update[:domain!] = Fᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Uii"
                update[:domain!] = Uᵢᵢ(CPU())
            end
        elseif instr[:fwrk][:deform] == "infinitesimal"
            if instr[:basis][:how] == "detΔFij"
                update[:domain!] = detΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔFii"
                update[:domain!] = ΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔUii"
                update[:domain!] = ΔUᵢᵢ(CPU())
            end
        end
    end
    if instr[:fwrk][:locking]
        update[:ΔJn!] = ΔJn(CPU())
        update[:ΔJs!] = ΔJs(CPU())
        update[:ΔJp!] = ΔJp(CPU())
    end
    
    if instr[:perf][:status]
        update[:elast!] = elast_fast(CPU())
    else
        update[:elast!] = elast(CPU())
    end


    update[:nonloc!] = nonlocal(CPU())
    if instr[:plast][:constitutive] == "MC"
        #ηmax = MCRetMap!(mpts,ϵpII,cmp,instr[:fwrk])
    elseif instr[:plast][:constitutive] == "DP"        
        if instr[:fwrk][:deform] == "finite"
            update[:retmap!] = finite_DP(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            update[:retmap!] = infinitesimal_DP(CPU())
        end
    elseif instr[:plast][:constitutive] == "J2"
        if instr[:fwrk][:deform] == "finite"
            update[:retmap!] = finite_J2(CPU())
        elseif instr[:fwrk][:deform] == "infinitesimal"
            update[:retmap!] = infinitesimal_J2(CPU())
        end
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mpts,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(instr[:plast][:constitutive])"))
    end
    
    return (;update...)
end

function elasto(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,dt::T2,instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
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

function elastoplast(mpts::Point{T1,T2,E,R},mesh::Mesh{T1,T2,D},cmpr::NamedTuple,dt::T2,instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
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

function thermo(mpts::Point{T1,T2,E,R},mesh::MeshThermalPhase{T1,T2},instr::Instruction{T1,T2,D}) where {T1,T2,E,R,D}
    # update temperature
    instr.cairn.update.heat!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    return nothing
end