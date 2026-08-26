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


    update[:nonloc_pq!] = nonlocal_pq(CPU())
    update[:nonloc_qp!] = nonlocal_qp(CPU())
    # dispatch resolves at kernel launch from Point's CM (DruckerPrager/VonMises) and ST
    # (LogarithmicStrain/InfinitesimalStrain) type parameters — see DP.jl's `retmap`
    # docstring. Unrecognized `plast.constitutive` strings now fail fast in `setup_cmp`
    # (setup time), not here.
    update[:retmap!] = retmap(CPU())

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
        for p ∈ 1:mpts.nmp
            mpts.s.ϵpII[p] = SVector{2,T2}(mpts.s.ϵpII[p][1], T2(0.0))
        end
        # CSR element→particle bucket list (O(nmp+nel), sequential) bounds nonlocal's
        # neighbor search to O(nmp×k) via basis.e2e, instead of the old O(nmp²) dense
        # e2p/p2p scan — see nonlocal.jl's docstring.
        ptr,idx = build_el2p(basis.p2e, T1(mesh.prprt.nel[end]))
        W       = zeros(T2,mpts.nmp)
        solver.cairn.update.nonloc_pq!(W,mpts,basis,T2(solver.nonloc.ls),ptr,idx; ndrange=mpts.nmp);sync(CPU())
        solver.cairn.update.nonloc_qp!(W,mpts,basis,T2(solver.nonloc.ls),ptr,idx; ndrange=mpts.nmp);sync(CPU())
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

