function elasto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::NamedTuple) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    instr[:cairn][:update].deform!(mpts,mesh.s,dt; ndrange=mpts.nmp);sync(CPU())
    # update material point's domain
    if instr[:basis][:which] == "gimpm"
        instr[:cairn][:update].domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    if instr[:fwrk][:locking]
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0)/mesh.prprt.dim
        # mapping to mesh 
        instr[:cairn][:update].ΔJn!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr[:cairn][:update].ΔJs!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
        # compute determinant Jbar 
        instr[:cairn][:update].ΔJp!(mpts,mesh,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    instr[:cairn][:update].elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())
    return nothing
end

function elastoplast(mpts::Point{T1,T2},mesh::Mesh{T1,T2},cmpr::NamedTuple,dt::T2,instr::NamedTuple) where {T1,T2}
    # update {logarithmic|infinitesimal} strains
    instr[:cairn][:update].deform!(mpts,mesh.s,dt; ndrange=mpts.nmp);sync(CPU())
    # update material point's domain
    if instr[:basis][:which] == "gimpm"
        instr[:cairn][:update].domain!(mpts; ndrange=mpts.nmp);sync(CPU())
    end
    # volumetric locking correction
    if instr[:fwrk][:locking]
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0)/mesh.prprt.dim
        # mapping to mesh 
        instr[:cairn][:update].ΔJn!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
        # compute nodal determinant of incremental deformation 
        instr[:cairn][:update].ΔJs!(mesh; ndrange=mesh.prprt.nno[end]);sync(CPU())
        # compute determinant Jbar 
        instr[:cairn][:update].ΔJp!(mpts,mesh,dim; ndrange=mpts.nmp);sync(CPU())
    end
    # update {kirchoff|cauchy} stresses
    instr[:cairn][:update].elast!(mpts,cmpr.Del; ndrange=mpts.nmp);sync(CPU())
    # plastic corrector
    if instr[:nonloc][:status]
        fill!(mpts.e2p,T1(0))
        fill!(mpts.p2p,T1(0))
        mpts.s.ϵpII[2,:].= T2(0.0)
        W,w     = spzeros(T2,mpts.nmp),spzeros(T2,mpts.nmp,mpts.nmp)
        for proc ∈ ["tplgy","p->q","p<-q"]
            instr[:cairn][:update].nonloc!(W,w,mpts,mesh,T2(instr[:nonloc][:ls]),proc; ndrange=mpts.nmp);sync(CPU())
        end
    else
        mpts.s.ϵpII[2,:].= mpts.s.ϵpII[1,:]
    end
    # plastic return-mapping dispatcher
    instr[:cairn][:update].retmap!(mpts,cmpr; ndrange=mpts.nmp);sync(CPU())
    return nothing
end

function thermo(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},instr::NamedTuple) where {T1,T2}
    # update temperature
    instr[:cairn][:update].heat!(mpts,mesh; ndrange=mpts.nmp);sync(CPU())
    return nothing
end