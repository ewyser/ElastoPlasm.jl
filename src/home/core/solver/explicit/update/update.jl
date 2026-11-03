function init_update(instr::NamedTuple; update::Dict{Symbol,Cairn} = Dict{Symbol,Cairn}())
    backend = exec_backend(instr)   # the solver's device (backend.exec), CPU() for the host
    update[:deform!] = deform(backend)
    
    update[:heat!] = heatflux(backend)
    if instr[:stab][:locking]
        update[:ΔJn!] = ΔJn(backend)
        update[:ΔJp!] = ΔJp(backend)
    end
    
    update[:elast!] = elast(backend)


    update[:nonloc_pq!] = nonlocal_pq(backend)
    update[:nonloc_qp!] = nonlocal_qp(backend)
    # the particles of each element, for the nonlocal neighbour search (element_particles, nonlocal.jl)
    update[:el2p_count!] = el2p_count(backend)
    update[:el2p_fill!]  = el2p_fill(backend)
    update[:el2p_sort!]  = el2p_sort(backend)
    # dispatch resolves at kernel launch from Point's CM (DruckerPrager/VonMises) and ST
    # (LogarithmicStrain/InfinitesimalStrain) type parameters — see DP.jl's `retmap`
    # docstring. Unrecognized `material.plastic` strings now fail fast in `setup_cmp`
    # (setup time), not here.
    update[:retmap!] = retmap(backend)

    # update mp's coordinates
    update[:coords!] = coord(backend)
    
    return (;update...)
end

function elasto(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},dt::T2,solver::Solver{T1,T2}) where {T1,T2}
    bk = exec_backend(solver)
    # update {logarithmic|infinitesimal} strains
    solver.cairn.update.deform!(mpts, mesh, basis, dt; ndrange=mpts.nmp)
    sync(bk)
    # volumetric locking correction
    if solver.stab.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0/(length(mesh.prprt.nel)-1))
        # mapping to mesh
        solver.cairn.update.ΔJn!(mpts, mesh, basis; ndrange=mpts.nmp)
        sync(bk)
        # compute determinant Jbar
        solver.cairn.update.ΔJp!(mpts, mesh, basis, dim; ndrange=mpts.nmp)
        sync(bk)
    end
    # update {kirchoff|cauchy} stresses
    solver.cairn.update.elast!(mpts; ndrange=mpts.nmp)
    sync(bk)
    # update mp's coordinates
    solver.cairn.update.coords!(mpts, mesh; ndrange=mpts.nmp)
    sync(bk)
    return nothing
end

function elastoplast(mpts::Point{T1,T2},mesh::Mesh{T1,T2},basis::Basis{T1,T2},dt::T2,solver::Solver{T1,T2}) where {T1,T2}
    bk = exec_backend(solver)
    # update {logarithmic|infinitesimal} strains
    solver.cairn.update.deform!(mpts, mesh, basis, dt; ndrange=mpts.nmp)
    sync(bk)
    # volumetric locking correction
    if solver.stab.locking
        # init mesh quantities to zero
        fill!(mesh.ΔJ,T2(0.0))
        # calculate dimensional cst.
        dim = T2(1.0/(length(mesh.prprt.nel)-1))
        # mapping to mesh
        solver.cairn.update.ΔJn!(mpts, mesh, basis; ndrange=mpts.nmp)
        sync(bk)
        # compute determinant Jbar
        solver.cairn.update.ΔJp!(mpts, mesh, basis, dim; ndrange=mpts.nmp)
        sync(bk)
    end
    # update {kirchoff|cauchy} stresses
    solver.cairn.update.elast!(mpts; ndrange=mpts.nmp)
    sync(bk)
    # plastic corrector
    if solver.nonloc.status
        map!(e -> SVector{2,T2}(e[1], T2(0.0)), mpts.s.ϵpII, mpts.s.ϵpII)   # element-wise: any device
        # the particles of each element bound nonlocal's neighbour search to O(nmp×k) via basis.e2e, instead of
        # an O(nmp²) scan (see nonlocal.jl's docstring)
        el2p = element_particles(mpts,mesh,basis,solver)
        W    = fill!(similar(mpts.Ω), T2(0.0))
        solver.cairn.update.nonloc_pq!(W, mpts, basis, T2(solver.nonloc.ls), el2p; ndrange=mpts.nmp)
        sync(bk)
        solver.cairn.update.nonloc_qp!(W, mpts, basis, T2(solver.nonloc.ls), el2p; ndrange=mpts.nmp)
        sync(bk)
    else
        map!(e -> SVector{2,T2}(e[1], e[1]), mpts.s.ϵpII, mpts.s.ϵpII)
    end
    # plastic return-mapping dispatcher
    solver.cairn.update.retmap!(mpts; ndrange=mpts.nmp)
    sync(bk)
    # update mp's coordinates
    solver.cairn.update.coords!(mpts, mesh; ndrange=mpts.nmp)
    sync(bk)
    return nothing
end

function thermo(mpts::Point{T1,T2},mesh::MeshThermalPhase{T1,T2},basis::Basis{T1,T2},solver::Solver{T1,T2}) where {T1,T2}
    bk = exec_backend(solver)
    # update temperature
    solver.cairn.update.heat!(mpts, mesh, basis; ndrange=mpts.nmp)
    sync(bk)
    return nothing
end

