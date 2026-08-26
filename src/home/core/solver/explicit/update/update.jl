function init_update(instr::NamedTuple; update::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    if instr[:perf][:status]
        update[:deform!] = deform_fast(CPU())
    else
        update[:deform!] = deform(CPU())
    end
    
    update[:heat!] = heatflux(CPU())
    if instr[:basis][:which] == "gimpm"
        #update[:domain!] = undeformed(CPU())
        update[:domain!] = Uᵢᵢ(CPU())
        #=
        update[:domain!] = undeformed(CPU())
        if instr[:strain][:deform] == "finite"
            if instr[:basis][:how] == "detFij"
                update[:domain!] = detFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Fii"
                update[:domain!] = Fᵢᵢ(CPU())
            elseif instr[:basis][:how] == "Uii#
                update[:domain!] = Uᵢᵢ(CPU())
            end
        elseif instr[:strain][:deform] == "infinitesimal"
            if instr[:basis][:how] == "detΔFij"
                update[:domain!] = detΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔFii"
                update[:domain!] = ΔFᵢᵢ(CPU())
            elseif instr[:basis][:how] == "ΔUii"
                update[:domain!] = ΔUᵢᵢ(CPU())
            end
        end
        =#
    end
    if instr[:stab][:locking]
        update[:ΔJn!] = ΔJn(CPU())
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
        if instr[:strain][:deform] == "finite"
            update[:retmap!] = finite_DP(CPU())
        elseif instr[:strain][:deform] == "infinitesimal"
            update[:retmap!] = infinitesimal_DP(CPU())
        end
    elseif instr[:plast][:constitutive] == "J2"
        if instr[:strain][:deform] == "finite"
            update[:retmap!] = finite_J2(CPU())
        elseif instr[:strain][:deform] == "infinitesimal"
            update[:retmap!] = infinitesimal_J2(CPU())
        end
    elseif instr[:plast][:constitutive] == "camC"
        #ηmax = camCRetMap!(mpts,cmp,instr[:fwrk])
    else
        throw(error("InvalidReturnMapping: $(instr[:plast][:constitutive])"))
    end

    # update mp's coordinates
    update[:coords!] = coord(CPU())
    
    return (;update...)
end

function elasto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},dt::T2,solver::ExplicitSolver{T1,T2}) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    solver.cairn.update.deform!(mpts,mesh,basis,dt; ndrange=mpts.nmp);sync(CPU())
    # volumetric locking correction
    if solver.stab.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0/(length(mesh.prprt.nel)-1))
        # mapping to mesh
        solver.cairn.update.ΔJn!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        # compute determinant Jbar
        solver.cairn.update.ΔJp!(mpts,mesh,basis,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update material point's domain
    if solver.basis.which == "gimpm"
        solver.cairn.update.domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    solver.cairn.update.elast!(mpts; ndrange=mpts.nmp);sync(CPU())
    # update mp's coordinates
    solver.cairn.update.coords!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    return nothing
end

function elastoplast(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},dt::T2,solver::ExplicitSolver{T1,T2}) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    solver.cairn.update.deform!(mpts,mesh,basis,dt; ndrange=mpts.nmp);sync(CPU())
    # update material point's domain
    if solver.basis.which == "gimpm"
        solver.cairn.update.domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    if solver.stab.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0/(length(mesh.prprt.nel)-1))
        # mapping to mesh
        solver.cairn.update.ΔJn!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
        # compute determinant Jbar
        solver.cairn.update.ΔJp!(mpts,mesh,basis,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    solver.cairn.update.elast!(mpts; ndrange=mpts.nmp);sync(CPU())
    # plastic corrector
    if solver.nonloc.status
        fill!(basis.e2p,T1(0))
        fill!(basis.p2p,T1(0))
        for p ∈ 1:mpts.nmp
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], T2(0.0))
        end
        W,w     = spzeros(T2,mpts.nmp),spzeros(T2,mpts.nmp,mpts.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            solver.cairn.update.nonloc!(W,w,mpts,mesh,basis,T2(solver.nonloc.ls),proc; ndrange=mpts.nmp);sync(CPU())
        end
    else
        for p ∈ 1:mpts.nmp
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], mpts.s.ϵpII[p][1])
        end
    end
    # plastic return-mapping dispatcher
    solver.cairn.update.retmap!(mpts; ndrange=mpts.nmp);sync(CPU())
    # update mp's coordinates
    solver.cairn.update.coords!(mpts,mesh; ndrange=(mpts.nmp));sync(CPU())
    return nothing
end

function thermo(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2},solver::ExplicitSolver{T1,T2}) where {T1,T2}
    # update temperature
    solver.cairn.update.heat!(mpts,mesh,basis; ndrange=mpts.nmp);sync(CPU())
    return nothing
end

